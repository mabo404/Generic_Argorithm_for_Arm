% main_pick_representatives_allgens.m
% 複数世代の CSV をまとめて読み込み → 代表個体を一括抽出
clear; clc;

% 世代CSVが置いてあるフォルダを選択
baseDir = uigetdir(pwd, 'generation_***_objs.csv があるフォルダを選択');
if isequal(baseDir,0), fprintf('キャンセルしました。\n'); return; end

% 対象ファイルを収集
files = dir(fullfile(baseDir, 'generation_*_objs.csv'));
if isempty(files)
    error('フォルダに generation_*_objs.csv が見つかりません。');
end

% すべて連結（どの世代の何行かも追跡）
F_all = []; gen_of_row = []; row_in_gen = []; file_of_row = strings(0,1);
for i = 1:numel(files)
    fn = fullfile(files(i).folder, files(i).name);
    F  = readmatrix(fn);
    if ~ismatrix(F), error('%s は行列ではありません。', files(i).name); end
    n = size(F,1);
    g = sscanf(files(i).name, 'generation_%d_objs.csv');  % 世代番号抽出
    F_all      = [F_all; F];
    gen_of_row = [gen_of_row; repmat(g, n, 1)];
    row_in_gen = [row_in_gen; (1:n)'];
    file_of_row= [file_of_row; repmat(string(files(i).name), n, 1)];
end

% 列名（必要に応じて合わせて下さい）
obj_names = ["J_tau","J_pos","J_ori","J_sc","J_smooth"];

% セレクタ（default_selectors.m を使用）
selectors = default_selectors();

% 代表抽出（可行=J_sc==0 でフィルタ）
[sel, T] = pick_representatives(F_all, ...
    'selectors', selectors, ...
    'obj_names', obj_names, ...
    'feasible_only', true, 'tol', 1e-12, ...
    'norm','percentile','plow',5,'phigh',95);

% 代表が属する世代・行・ファイル名を付与
T.gen        = gen_of_row(T.index);
T.row_in_gen = row_in_gen(T.index);
T.file       = file_of_row(T.index);
T = movevars(T, {'gen','row_in_gen','file'}, 'Before', 'index');

% 出力
outcsv = fullfile(baseDir, 'allgens_representatives.csv');
writetable(T, outcsv);
save(fullfile(baseDir, 'allgens_representatives.mat'), 'T','sel','selectors','obj_names');

disp(T);
fprintf('[OK] 代表一覧を保存: %s\n', outcsv);

% ワークスペースに配置（任意）
assignin('base','T_allgens',T);
assignin('base','sel_allgens',sel);
