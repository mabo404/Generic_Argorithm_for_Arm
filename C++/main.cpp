#include <pagmo/algorithms/nsga2.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/population.hpp>
#include <pagmo/utils/multi_objective.hpp>       // fast_non_dominated_sorting

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using pagmo::vector_double;
namespace fs = std::filesystem;

// ====== ユーザー設定 ==================================
static constexpr int    POP_SIZE   = 48;
static constexpr int    MAX_GEN    = 100;
static constexpr double INVALID_FITNESS = 1e5; // しきい値（>1e5を不正扱い）
static constexpr double L1_FIXED   = 0.10;
static constexpr double L_LOW[4]   = {0.10,0.10,0.10,0.10};      // L2..L5
static constexpr double L_UP [4]   = {0.30,0.30,0.20,0.20};      // L2..L5
// NSGA-II（MatlabのSBX/多項式変異パラメータに対応）
static constexpr double CX_RATE    = 0.99;     // 交叉確率（SBX）
static constexpr double ETA_C      = 1.0;     // SBX 分布指数
static constexpr double MUT_RATE   = 1.0/5.0; // 多項式変異 確率（各遺伝子）
static constexpr double ETA_M      = 2.0;     // 多項式変異 分布指数
static constexpr unsigned SEED     = 42u;     // 再現性

// ====== 便利関数 ============================================================
std::string timestamp_yyyymmdd_hhmmss() {
    using clock = std::chrono::system_clock;
    auto now = clock::now();
    std::time_t t = clock::to_time_t(now);
    std::tm tm{};
#ifdef _WIN32
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
    return oss.str();
}

struct Csv {
    std::ofstream ofs;
    explicit Csv(const fs::path &p, const std::string &header = "") {
        ofs.open(p, std::ios::out);
        if (!ofs) throw std::runtime_error("Cannot open: " + p.string());
        if (!header.empty()) ofs << header << "\n";
    }
    void row(const std::vector<std::string> &cells) {
        bool first = true;
        for (auto &c : cells) {
            if (!first) ofs << ",";
            ofs << c;
            first = false;
        }
        ofs << "\
        n";
    }
};

// ====== あなたの「評価関数」を入れる場所（TODO: 後でEigen版に差し替え） =====
static inline vector_double evaluate_arm_dummy(const vector_double &x) {
    // x = [L1..L5]
    // --- ダミー: 形だけ4目的を返す（すべて有限値にする）
    // 実装方針: 後でここを Eigen を使った本物の Evaluate(...) に置き換える
    const double Lsum = x[0]+x[1]+x[2]+x[3]+x[4];
    const double L2   = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4];
    const double J_tau =  0.5*Lsum;        // 例: 小さいほど良い
    const double J_pos =  L2;              // 例
    const double J_ori =  (x[1]+x[3]);     // 例
    const double J_sc  =  (x[2]+x[4]);     // 例
    return {J_tau, J_pos, J_ori, J_sc};
}

// ====== UDP: pagmo に渡す「問題クラス」 ====================================
struct ArmUDP {
    // 目的関数（4目的を返す）
    pagmo::vector_double::size_type get_nobj() const { return 4u; }
    vector_double fitness(const vector_double &x) const {
        // x は L2..L5 の4次元。L1を内部で固定して5次元ベクトルに再構成
        vector_double full_x = {L1_FIXED, x[0], x[1], x[2], x[3]};
        auto f = evaluate_arm_dummy(full_x);
        for (auto &v : f) {
            if (!std::isfinite(v) || v > INVALID_FITNESS) v = INVALID_FITNESS;
        }
        return f;
    }
    // 変数の下限・上限
    std::pair<vector_double, vector_double> get_bounds() const {
        vector_double lb(4), ub(4);
        for (int i=0;i<4;++i){ lb[i]=L_LOW[i]; ub[i]=L_UP[i]; }
        return {lb, ub};
    }
    // バッチ評価（並列化したくなったら有効化）
    // bool has_batch_fitness() const { return true; }
    // std::vector<vector_double> batch_fitness(const std::vector<vector_double> &xs) const {
    //     std::vector<vector_double> out(xs.size());
    //     #pragma omp parallel for
    //     for (int i=0;i<(int)xs.size();++i) out[i] = fitness(xs[i]);
    //     return out;
    // }
};

// ====== 指標の集計（Matlabのログに対応） ===================================
struct GenStats {
    double gen_time_sec{};
    int invalid_count{};
    int viol_safety{}; // J_tau > 0 の個体数
    vector_double best_obj; // ここでは J_tau を最小にする個体の目的ベクトル
    vector_double best_vars; // その設計変数
};

// 集団から統計を作る
GenStats make_stats(const pagmo::population &pop, double elapsed_sec) {
    GenStats s;
    s.gen_time_sec = elapsed_sec;
    const auto &F = pop.get_f();
    const auto &X = pop.get_x();
    int n = static_cast<int>(F.size());
    s.invalid_count = 0;
    s.viol_safety   = 0;
    int best_idx = 0;
    double best_Jtau = std::numeric_limits<double>::infinity();
    for (int i=0;i<n;++i) {
        bool valid = true;
        for (double v : F[i]) if (!std::isfinite(v) || v > INVALID_FITNESS) { valid=false; break; }
        if (!valid) s.invalid_count++;
        if (valid && F[i][0] > 0.0) s.viol_safety++; // J_tau>0
        if (F[i][0] < best_Jtau) { best_Jtau = F[i][0]; best_idx = i; }
    }
    s.best_obj  = F[best_idx];
    s.best_vars = { L1_FIXED,                   // ここで5次元に拡張
        X[best_idx][0], X[best_idx][1], X[best_idx][2], X[best_idx][3] };
    return s;
}

int main() try {
    // 出力フォルダ run_YYYYMMDD_HHMMSS/
    const auto ts = timestamp_yyyymmdd_hhmmss();
    const fs::path logdir = fs::current_path() / ("run_" + ts);
    fs::create_directories(logdir);

    // CSV: デバッグサマリ（Matlabのdebug_log.csv相当）
    Csv log(logdir / "debug_log.csv",
            "gen,gen_time_sec,invalid_count,viol_safety,"
            "best_Jtau,best_Jpos,best_Jori,best_Jsc,"
            "best_L1,best_L2,best_L3,best_L4,best_L5");

    // 問題＆アルゴリズム＆初期集団
    pagmo::problem prob{ArmUDP{}};
    // gen=1 にして自前ループでログを毎世代出せるようにする
    pagmo::algorithm algo{pagmo::nsga2(1, CX_RATE, ETA_C, MUT_RATE, ETA_M, SEED)};
    algo.set_seed(SEED);
    pagmo::population pop{prob, POP_SIZE, SEED};

    // 初期評価（ベースライン: gen=0）
    auto t0 = std::chrono::high_resolution_clock::now();
    auto t1 = std::chrono::high_resolution_clock::now();
    GenStats s0 = make_stats(pop, std::chrono::duration<double>(t1-t0).count());
    log.row(std::vector<std::string>{ "0",
              std::to_string(s0.gen_time_sec),
              std::to_string(s0.invalid_count),
              std::to_string(s0.viol_safety),
              std::to_string(s0.best_obj[0]), std::to_string(s0.best_obj[1]),
              std::to_string(s0.best_obj[2]), std::to_string(s0.best_obj[3]),
              std::to_string(s0.best_vars[0]), std::to_string(s0.best_vars[1]),
              std::to_string(s0.best_vars[2]), std::to_string(s0.best_vars[3]),
              std::to_string(s0.best_vars[4]) });

    // 世代ループ
    double first_gen_wall = 0.0;
    for (int gen=1; gen<=MAX_GEN; ++gen) {
        auto tg0 = std::chrono::high_resolution_clock::now();
        pop = algo.evolve(pop); // 1世代だけ進化
        auto tg1 = std::chrono::high_resolution_clock::now();
        double dt = std::chrono::duration<double>(tg1 - tg0).count();

        // 1世代目で総見積りを表示（Matlabに合わせて）
        if (gen == 1) {
            first_gen_wall = dt;
            const double eval_time_per_ind = first_gen_wall / (2.0 * POP_SIZE);
            const double total_eval_count  = (MAX_GEN - 1) * POP_SIZE;
            const double est_total_sec     = eval_time_per_ind * total_eval_count;
            int hrs = int(est_total_sec / 3600);
            int mins = int(std::fmod(est_total_sec,3600.0) / 60.0);
            int secs = int(std::fmod(est_total_sec,60.0));
            std::cout << "第1世代にかかった時間 = " << first_gen_wall
                      << " 秒（" << 2*POP_SIZE << "体分）\n";
            std::cout << "全体の予測実行時間（評価体数に基づく） = 約 "
                      << hrs << " 時間 " << mins << " 分 " << secs
                      << " 秒（" << est_total_sec << " 秒）\n";
        }

        // 集計＆ログ
        GenStats sg = make_stats(pop, dt);
        log.row(std::vector<std::string>{ std::to_string(gen),
                  std::to_string(sg.gen_time_sec),
                  std::to_string(sg.invalid_count),
                  std::to_string(sg.viol_safety),
                  std::to_string(sg.best_obj[0]), std::to_string(sg.best_obj[1]),
                  std::to_string(sg.best_obj[2]), std::to_string(sg.best_obj[3]),
                  std::to_string(sg.best_vars[0]), std::to_string(sg.best_vars[1]),
                  std::to_string(sg.best_vars[2]), std::to_string(sg.best_vars[3]),
                  std::to_string(sg.best_vars[4]) });
    }

    // 最終：第1フロントを抽出してCSVに保存
    const auto &F = pop.get_f();
    const auto &X = pop.get_x();
    auto [fronts, ranks, dom_count, dominated_by] = pagmo::fast_non_dominated_sorting(F); // fronts[0] が第1フロント
    Csv pareto_objs (logdir / "ParetoObjs.csv", "J_tau,J_pos,J_ori,J_sc");
    Csv pareto_vars (logdir / "ParetoVars.csv", "L1,L2,L3,L4,L5");
    for (auto idx : fronts[0]) {
        pareto_objs.row(std::vector<std::string>{ std::to_string(F[idx][0]), std::to_string(F[idx][1]),
                          std::to_string(F[idx][2]), std::to_string(F[idx][3]) });
        const auto &v = X[idx]; // v は4次元 (L2..L5)
        pareto_vars.row(std::vector<std::string>{
            std::to_string(L1_FIXED),
            std::to_string(v[0]), std::to_string(v[1]), std::to_string(v[2]), std::to_string(v[3])
        });
    }

    std::cout << "done. results in: " << logdir << "\n";
    return 0;

} catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
}
