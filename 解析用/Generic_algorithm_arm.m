% === 概要 ===============================================================
% 目的 : 6軸アームのリンク長の多目的GA最適化（NSGA-II系）
% 入力 : なし（本スクリプト内で乱数・境界・世代数を定義）
% 出力 : run_YYYYMMDD_HHMMSS/ に CSV / MAT を保存
% 仕様 : 目的関数 [J_tau, J_pos, J_ori, J_sc]（列順固定）
% 単位 : 角度[rad], 長さ[m], トルク[N·m]（座標系はW=世界, B=ベース）
% 再現 : rng(seed_cpu, alg_cpu) を固定
% =======================================================================

%% 01. 初期設定とメタ情報
% - 乱数・日付・出力ディレクトリ名など
clearvars; % ワークスペース内の変数を全消去（関数・パス・ブレークポイントには影響なし）
clc; % コマンドウィンドウの表示を消去（履歴のみクリア、変数は保持）
close all; % 開いているすべての Figure ウィンドウを閉じる

%% 02. ハイパーパラメータと設計変数の境界
% - GA設定（個体数/世代数/トーナメントk など）
% - SBX/多項式変異のパラメータ
% - 設計変数の下限/上限
INVALID_FITNESS = 1e5; % 評価値の無効判定しきい値（>1e5 を不正扱い。NaN/Inf 判定と併用）
num_vari = 5; % 設計変数の数（L1〜L5の計5）
dis_c = 1; % SBX（交叉）の分布指数（小さいほど親に近い子が出やすい）
dis_m = 2; % 多項式変異の分布指数（大きいほど変異幅が小さく局所的になる）
pro_m = 1/num_vari; % 各遺伝子の変異確率（1/変数数）—平均で1個体あたり約1遺伝子が変異
pop_size = 2; % 集団サイズ（1世代あたりの個体数）※SBXの親ペア生成のため偶数推奨
assert(mod(pop_size,2) == 0, 'SBX のため、pop_size は偶数でなければなりません。'); % SBXの親ペア形成に必要なため、集団サイズが偶数であることを検証（奇数ならエラー）
max_gen = 1; % 進化させる総世代数（ループ回数）
minL = 0.1; % 各リンクの最小長さ [m]（L1〜L5 の下限）
lower_bound = [minL, minL, minL, minL, minL]; % 設計変数の下限 [L1..L5]（Lは minL）[m]
upper_bound = [0.1,  0.3,  0.3,  0.2,  0.2]; % 設計変数の上限 [L1..L5]（Lは0.5[m]）
firstGen_start = tic; % 第1世代の処理時間計測を開始（タイマ開始）
generation = 1; % 現在の世代番号を初期化（1から開始）
k = 2; % トーナメント選択の参加者数 k（各トーナメントでの比較個体数）

%% 03. ログの準備
% - ログ用フォルダ作成、CSV/デバッグファイルのパスを決定
% - （必要なら）ヘッダ行の初期化
ts = string(datetime('now','Format','yyyyMMdd_HHmmss')); % 実行時刻のタイムスタンプ文字列を生成（例: 20250820_142755）—ログ/出力ディレクトリ名に使用
logdir = char(fullfile(pwd, "run_" + ts)); % カレントフォルダ直下の出力先パス（'run_タイムスタンプ'）を組み立てる※まだ作成はしない
if ~isfolder(logdir) % 出力用ディレクトリが存在しなければ作成（初回実行時の安全策）
    mkdir(logdir);
end
debug_log = fullfile(logdir,'debug_log.csv'); % 書き込みは log_write_generation が行う
seed_cpu = 42; % 乱数生成の種（CPU側）。再現性確保のため固定値にする
alg_cpu = 'twister'; % 乱数生成アルゴリズム（'twister'＝Mersenne Twister）を指定して再現性を統一
rng(seed_cpu, alg_cpu); % 乱数ストリームを初期化（種とアルゴリズムで再現性を固定）

%% 04. 初期集団の生成
% - 境界内で初期個体を生成（必要なら既知解を混ぜる）
pop_vari = repmat(lower_bound,pop_size,1) ...
           + rand(pop_size,num_vari) .* repmat((upper_bound-lower_bound),pop_size,1); % 下限＋一様乱数×(上限−下限)で境界内の初期個体を生成（行=個体, 列=設計変数）

%% 05. 初期評価とベースライン記録
% - 初期集団を評価し、フロント抽出とベスト個体を記録
pop_fitness = Evaluate(pop_vari); % 初期集団の目的関数を評価（行=個体，列=[J_tau, J_pos, J_ori, J_sc]）
num_obj = size(pop_fitness, 2); % 目的関数の数（列数）を取得—以降の配列サイズ確保に使用
gen_time_sec = zeros(max_gen,1); % 各世代の処理時間 [s] を格納するベクトルを事前確保（max_gen×1）
invalid_count = zeros(max_gen,1); % 各世代の不正個体数（無限大/NaN/しきい値超え）を記録する配列を初期化
viol_safety = zeros(max_gen,1); % 各世代の安全違反数を記録する配列を初期化（例: J_tau>0 を違反とみなす）
rank_counts = zeros(max_gen,10); % 各世代のランク別個体数を記録（列=ランク）。まず10列ぶん確保し、必要なら後で列を拡張。
best_obj_record = zeros(max_gen+1, num_obj); % ベスト目的値の履歴（初期ベースライン用に+1行確保）。行=世代0..max_gen、列=各目的関数
param_history = zeros(max_gen+1, num_vari); % ベスト個体の設計変数履歴（初期ベースライン分+1行）。行=世代0..max_gen, 列=設計変数
rank_now = GA.pareto_rank(pop_fitness); % 初期集団に対するパレートランクを計算（1が最良フロント）
front_idx = find(rank_now == 1); % 第1フロント（ランク=1）の個体インデックスを抽出
F_front = pop_fitness(front_idx, :); % 第1フロントの目的関数値を抽出（行=該当個体, 列=[J_tau,J_pos,J_ori,J_sc]）
if isscalar(front_idx)
    idx0 = front_idx;
else
    cd0 = GA.crowding_distance(F_front);
    [~, rel0] = max(cd0);
    idx0 = front_idx(rel0);
end
best_obj = pop_fitness(idx0, :); % 代表個体 idx0 の目的関数ベクトルを取得（[J_tau, J_pos, J_ori, J_sc]）
best_vars = pop_vari(idx0, :); % 代表個体 idx0 の設計変数ベクトルを取得（[L1..L5]）
best_obj_record(1,:) = best_obj; % 初期ベースラインとして、世代0行（1行目）にベスト目的値を保存
param_history(1,:) = best_vars; % 初期ベースラインとして、世代0行（1行目）にベスト個体の設計変数を保存

%% 06. 世代ループ（GA本体）
% - 選択 → 交叉 → 変異 → 子評価 → 環境選択 → ログ出力
while generation < max_gen % 世代ループの開始条件（現在世代が最大世代に達するまで反復）
    tgen = tic; % この世代の処理時間の計測を開始（タイマ開始）

    % ===== [A] 親個体の選択（トーナメント + クラウディング距離） =====
    tmp = randi(pop_size, pop_size, k); % 各行にk個の候補個体IDを持つ行列を生成（1..pop_size から乱数選択）
    rank = GA.pareto_rank(pop_fitness); % 現在の集団に対してパレートランクを計算（1が最良フロント）
    cd_all = GA.compute_crowding_all(pop_fitness, rank); % フロントごとに算出したクラウディング距離を、全個体分のベクトルとして集約する
    sel = GA.tournament_select(tmp, rank, cd_all); % トーナメント選択で親個体のインデックスを決定（優先: ランク最小、同率はクラウディング距離が大きい方）
    parents = pop_vari(sel,:); % 選択された親個体の設計変数行列を抽出（行=親, 列=設計変数）
    
    % ===== [B] 交叉（SBX） =====
    pop_crossover = GA.sbx_crossover(parents, lower_bound, upper_bound, dis_c); % SBX（模擬二進交叉）で子個体を生成（境界を尊重、分布指数 dis_c を使用）

    % ===== [C] 変異（多項式変異） =====
    pop_mutation = GA.poly_mutation(pop_crossover, lower_bound, upper_bound, dis_m, pro_m); % 多項式変異で子個体に摂動を付与（分布指数 dis_m・変異確率 pro_m を使用）

    % ===== [D] 子個体の評価 =====
    pop_mutation_fitness = Evaluate(pop_mutation); % 生成した子個体の目的関数を評価（行=個体，列=[J_tau,J_pos,J_ori,J_sc]）

    % ===== [E] 環境選択（親+子 → 次世代へ選抜） =====
    combined_vari = [pop_vari; pop_mutation]; % 親と子の設計変数を縦に結合（次段の環境選択に備える）
    combined_fit = [pop_fitness; pop_mutation_fitness]; % 親と子の目的関数を縦に結合（環境選択の入力とする）
    [pop_vari, pop_fitness] = GA.environmental_selection(combined_vari, combined_fit, pop_size); % 環境選択で合体集団から次世代を選抜（NSGA-II：ランク優先→クラウディング距離）
    gen_time_sec(generation) = toc(tgen); % この世代の経過時間 [s] を記録（タイマ停止）
    invalid_count(generation) = sum(any(pop_fitness > INVALID_FITNESS | ~isfinite(pop_fitness), 2)); % 無効個体数をカウント（しきい値超え or NaN/Inf がある行を合計）
    viol_safety(generation) = sum(pop_fitness(:,1) > 0); % 安全違反数を集計（J_tau>0 を違反と定義し、第1列が正の個体をカウント）
    rk = GA.pareto_rank(pop_fitness); % 次世代に対するパレートランクを再計算（1が最良）—統計や代表抽出に使用
    K = max(rk); % 最大ランク値＝フロント数を取得（1..K）
    if size(rank_counts,2) < K % ランク数が既存の列数を超える場合は列拡張が必要か判定
        rank_counts(:, end+1:K) = 0; % 足りない列（end+1〜K）を0で埋めて確保
    end
    rank_counts(generation, 1:K) = accumarray(rk, 1, [K,1]).'; % 当該世代のランク別個体数を集計して1..K列に格納（accumarrayで各ランクをカウントし、転置して行ベクトル化）
    front_idx = find(rk == 1); % 第1フロントのインデックスを抽出（ランク=1の個体集合）
    F_front = pop_fitness(front_idx, :); % 第1フロント個体の目的関数行列を抽出（行=該当個体, 列=[J_tau,J_pos,J_ori,J_sc]）
    if isscalar(front_idx) % 第1フロントが1個体なら距離比較せずその個体を採用
        bestIdx = front_idx; % ベスト個体のインデックスにその単一個体を設定
    else % （それ以外の場合）複数個体から代表を選ぶ
        cd_front = GA.crowding_distance(F_front); % 第1フロント内のクラウディング距離を計算
        [~, bestRelIdx] = max(cd_front); % 距離が最大の相対インデックスを取得（多様性最大を優先）
        bestIdx = front_idx(bestRelIdx); % 絶対インデックスに変換してベスト個体を確定
    end
    best_obj = pop_fitness(bestIdx,:); % ベスト個体の目的関数ベクトルを取得（[J_tau, J_pos, J_ori, J_sc]）
    best_vars = pop_vari(bestIdx,:); % ベスト個体の設計変数ベクトルを取得（[L1..L5]）
    best_obj_record(generation+1,:) = best_obj; % ベスト目的値を履歴に保存（初期ベースラインが1行目のため，世代gはg+1行目に記録）
    param_history(generation+1,:) = best_vars; % ベスト個体の設計変数を履歴に保存（初期ベースライン分+1オフセットで格納）
    
    % ===== [F] 実行時間の予測 =====
    if generation == 1 % 第1世代終了時のみ、所要時間から総見積りを計算
        firstGen_time = toc(firstGen_start); % 第1世代終了までの経過時間 [s]（初期評価も含む）
        eval_time_per_individual = firstGen_time / (2 * pop_size); % 個体1体あたり時間の推定（※この時点で約2*pop_size体評価済みと仮定）
        total_eval_count = (max_gen - 1) * pop_size; % 残り世代で評価する総個体数（各世代で子 pop_size 体を想定）
        estimated_total_time = eval_time_per_individual * total_eval_count; % 残り総所要時間 [s] の見積り
        total_sec = round(estimated_total_time); % 見積り秒数を四捨五入
        hrs = floor(total_sec / 3600); % 時間成分
        mins = floor(mod(total_sec, 3600) / 60); % 分成分
        secs = mod(total_sec, 60); % 秒成分
        fprintf('第1世代にかかった時間 = %.2f 秒（%d体分）\n', firstGen_time, 2*pop_size); % 第1世代までの実測時間を表示
        fprintf('全体の予測実行時間（評価体数に基づく） = 約 =%d 時間 %d 分 %d 秒（%.2f 秒）\n', hrs, mins, secs, estimated_total_time); % 人が読みやすい形で残り時間を表示
    end

    % ===== [E] メタ計算（この世代） =====
    mask_valid  = all(isfinite(pop_fitness),2) & all(pop_fitness <= INVALID_FITNESS,2);
    invalid_count_gen = sum(~mask_valid);
    viol_safety_gen = sum(mask_valid & (pop_fitness(:,1) > 0));
    [best_obj_gen, best_idx_gen] = min(pop_fitness(:,1), [], 'omitnan');
    best_vars_gen = pop_vari(best_idx_gen, :);

    % ===== [F] ログ出力（CSV/デバッグ） =====
    meta_loop = struct( ...
        'gen_time_sec',  gen_time_sec(generation), ...
        'invalid_count', invalid_count(generation), ...
        'viol_safety',   viol_safety(generation), ...
        'best_obj',      best_obj_gen, ...       % ← 直前で計算済みを使用
        'best_vars',     best_vars_gen ...       % ← 同上
    );
    Log.log_write_generation(generation, pop_fitness, pop_vari, logdir, meta_loop);

    generation = generation + 1;   % ← while の最後でインクリメント
end

%% 07. 後処理と最終保存
% - 最終フロント抽出、メタ情報整理、results.mat の保存
rank_final = GA.pareto_rank(pop_fitness); % 最終世代のパレートランクを計算（1が最良フロント）
front1_idx = find(rank_final == 1); % 第1フロント（ランク=1）の個体インデックスを抽出
ParetoVars = pop_vari(front1_idx,:); % 第1フロント個体の設計変数行列を取得（行=個体, 列=設計変数）
ParetoObjs = pop_fitness(front1_idx,:); % 第1フロント個体の目的関数行列を取得（行=個体, 列=[J_tau,J_pos,J_ori,J_sc]）

meta = struct(); % 実行メタ情報を格納する構造体を初期化
meta.timestamp = ts; % 実行タイムスタンプ（run名と一致）
meta.logdir = logdir; % 出力ディレクトリのパス
meta.console = fullfile(logdir,'debug_log.csv'); % デバッグサマリCSV（console相当）のパス
meta.seed_cpu = seed_cpu; % 乱数シード
meta.alg_cpu = alg_cpu; % 乱数アルゴリズム名（'twister'等）
meta.pop_size = pop_size; % 集団サイズ
meta.max_gen = max_gen; % 最大世代数
meta.k_tour = k; % トーナメントの参加者数 k
meta.dis_c = dis_c; % SBX の分布指数
meta.dis_m = dis_m; % 多項式変異の分布指数
meta.pro_m = pro_m; % 変異確率
meta.lower_bound = lower_bound; % 設計変数の下限ベクトル [1x5]
meta.upper_bound = upper_bound; % 設計変数の上限ベクトル [1x5]
meta.matlab_ver = version; % MATLAB 本体のバージョン文字列
meta.toolbox_ver = ver; % インストール済みツールボックス情報（構造体配列）

try % 実行中のスクリプト/関数のフルパス取得を試みる
    meta.script = mfilename('fullpath'); % 取得できた場合はパスを格納
catch % 取得できない実行形態（無名スクリプト等）の場合
    meta.script = ''; % 空文字をセットしてフォールバック
end % try-catch 終了

save(fullfile(logdir,'results.mat'), 'meta', 'best_obj_record','param_history',...
'gen_time_sec','invalid_count','viol_safety','rank_counts', 'ParetoVars','ParetoObjs',...
'front1_idx', 'pop_vari','pop_fitness','rank_final'); % 出力フォルダに results.mat を作成し、メタ情報・履歴・最終集団・第1フロントなど主要変数一式を保存（再解析用）