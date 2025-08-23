% This is a MATLAB implementation of real-coded GA algorithm[1][2]
% [1] K. Deb, A. Kumar. Real-coded genetic algorithms 
%     with simulated binary crossover: studies on multimodal and multiobjective 
%     problems. Complex Systems, 1995, 9(6):431-54.
% [2] K. Deb. An efficient constraint handling method
%     for genetic algorithms. Computer Methods in Applied Mechanics and Engineering, 2000, 186(2):311-38.
%------------------------------------
clearvars;   % ワークスペースの変数をすべてクリアする
clc;         % コマンドウィンドウをクリアする
close all;   % 全ての図ウィンドウを閉じる
%------------------------------------
% 問題の選択: Ellipsoid, Rosenbrock, Ackley, Griewank のいずれか
obj_fun = 'Fun_Ellipsoid';   % 最適化対象の関数名を文字列で指定する
% 変数の次元数設定
num_vari = 50;               % デザイン変数の個数（次元）を50に設定する
% 集団サイズと世代数の設定
pop_size = 100;              % 各世代の個体数を100に設定する
max_gen  = 500;              % 最大世代数を500に設定する
%------------------------------------
% 各関数ごとの探索空間境界を定義
switch obj_fun
    case 'Fun_Ellipsoid'
        lower_bound = -5.12 * ones(1, num_vari);  % 各変数の下限を-5.12に設定
        upper_bound =  5.12 * ones(1, num_vari);  % 上限を+5.12に設定
    case 'Fun_Rosenbrock'
        lower_bound = -2.048 * ones(1, num_vari);
        upper_bound =  2.048 * ones(1, num_vari);
    case 'Fun_Ackley'
        lower_bound = -32.768 * ones(1, num_vari);
        upper_bound =  32.768 * ones(1, num_vari);
    case 'Fun_Griewank'
        lower_bound = -600 * ones(1, num_vari);
        upper_bound =  600 * ones(1, num_vari);
end
% 最良目的関数値を世代ごとに記録する配列を初期化
best_obj_record = zeros(max_gen, 1);
%------------------------------------
% 第1世代の生成と評価
generation = 1;  % 現在の世代を1に設定
% 一様乱数で初期集団を生成し、各変数を境界内にスケーリング
pop_vari   = repmat(lower_bound, pop_size, 1) + rand(pop_size, num_vari) .* repmat((upper_bound - lower_bound), pop_size, 1);
pop_fitness = feval(obj_fun, pop_vari);  % 各個体の目的関数値を計算
best_obj_record(generation) = min(pop_fitness);  % 第1世代の最良値を記録
% グラフ表示とコンソール出力
plot(best_obj_record(1:generation));  % 最良値の推移をプロット
xlabel('generation'); ylabel('best objective value');
title(sprintf('GA on %d-d %s function\n generation: %d, best: %0.4g', num_vari, obj_fun(5:end), generation, best_obj_record(generation))); drawnow;
fprintf('GA on %s, generation: %d, evaluation: %d, best: %0.4g\n', obj_fun(5:end), generation, generation * pop_size, best_obj_record(generation));
%------------------------------------
% 世代進化ループ
while generation < max_gen
    % 親の選択：2-トーナメント選択
    k     = 2;  % トーナメントサイズを2に設定
    temp  = randi(pop_size, pop_size, k);  % 各トーナメントでランダムに2個体を選ぶ
    [~, index] = min(pop_fitness(temp), [], 2);  % 目的関数が小さい方を勝者とする
    pop_parent = pop_vari(sub2ind(size(pop_vari), temp(:,1) .* (index == 1) + temp(:,2) .* (index == 2), ones(pop_size,1)), :);
    % 交叉：Simulated Binary Crossover (SBX)
    dis_c   = 1;  % SBX の分布指数
    mu       = rand(pop_size/2, num_vari);  % SBX 用乱数行列
    parent1  = pop_parent(1:2:end, :);  % 偶数番目の親群
    parent2  = pop_parent(2:2:end, :);  % 奇数番目の親群
    beta     = 1 + 2 * min(parent1, parent2) - lower_bound ./ max(abs(parent2 - parent1), 1e-6);
    alpha    = 2 - beta.^(-dis_c - 1);
    betaq    = (alpha .* mu).^(1/(dis_c + 1)) .* (mu <= 1./alpha) + (1./(2 - alpha .* mu)).^(1/(dis_c + 1)) .* (mu > 1./alpha);
    betaq    = betaq .* ((-1) .^ randi([0,1], pop_size/2, num_vari));  % SBX の符号付変更
    offspring1 = 0.5 * ((1 + betaq) .* parent1 + (1 - betaq) .* parent2);  % 子1
    offspring2 = 0.5 * ((1 - betaq) .* parent1 + (1 + betaq) .* parent2);  % 子2
    pop_crossover = [offspring1; offspring2];  % 生成子集団を結合
    % 突然変異：多項式突然変異
    dis_m    = 1;       % 変異の分布指数
    pro_m    = 1/num_vari;  % 各遺伝子の変異確率
    rand_var = rand(pop_size, num_vari);  % 変異判定用乱数
    mu       = rand(pop_size, num_vari);  % 変異量計算用乱数
    deta     = min(pop_crossover - lower_bound, upper_bound - pop_crossover) ./ (upper_bound - lower_bound);
    detaq    = zeros(pop_size, num_vari);
    pos1     = rand_var <= pro_m & mu <= 0.5;  % 低側変異条件
    pos2     = rand_var <= pro_m & mu >  0.5;  % 高側変異条件
    detaq(pos1) = (2 * mu(pos1) + (1 - 2 * mu(pos1)) .* (1 - deta(pos1)).^(dis_m + 1)).^(1/(dis_m + 1)) - 1;
    detaq(pos2) = 1 - (2 * (1 - mu(pos2)) + 2 * (mu(pos2) - 0.5) .* (1 - deta(pos2)).^(dis_m + 1)).^(1/(dis_m + 1));
    pop_mutation = pop_crossover + detaq .* (upper_bound - lower_bound);
    % 突然変異後の適応度計算
    pop_mutation_fitness = feval(obj_fun, pop_mutation);
    % 環境選択：親+子を評価順にソートし上位を次世代とする
    pop_combined_vari   = [pop_vari; pop_mutation];
    pop_combined_fitness = [pop_fitness; pop_mutation_fitness];
    [~, win_idx] = sort(pop_combined_fitness);
    pop_vari   = pop_combined_vari(win_idx(1:pop_size), :);
    pop_fitness = pop_combined_fitness(win_idx(1:pop_size));
    % 世代更新と記録
    generation = generation + 1;
    best_obj_record(generation) = min(pop_fitness);
    % 結果表示更新
    plot(best_obj_record(1:generation)); drawnow;
    fprintf('GA on %s, generation: %d, evaluation: %d, best: %0.4g\n', obj_fun(5:end), generation, generation * pop_size, best_obj_record(generation));
end
%------------------------------------
% 目的関数定義：Ellipsoid
function f = Fun_Ellipsoid(x)
    f = sum((1:size(x,2)) .* x.^2, 2);  % 重み付き2乗和
end
% Rosenbrock
function f = Fun_Rosenbrock(x)
    f = sum(100 * (x(:,2:end) - x(:,1:end-1).^2).^2 + (x(:,1:end-1) - 1).^2, 2);
end
% Ackley
function f = Fun_Ackley(x)
    d = size(x,2);
    f = -20 * exp(-0.2 * sqrt(sum(x.^2,2)/d)) - exp(sum(cos(2*pi*x),2)/d) + 20 + exp(1);
end
% Griewank
function f = Fun_Griewank(x)
    d = size(x,2);
    f = sum(x.^2/4000,2) - prod(cos(x./sqrt(1:d)),2) + 1;
end