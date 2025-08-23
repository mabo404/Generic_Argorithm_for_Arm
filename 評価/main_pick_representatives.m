% main_pick_representatives.m
% 目的行列(N×M)ファイル（CSV/MAT）を選択 → 代表個体を抽出 → 表示＆保存
clear; clc;

[fn, fp] = uigetfile({'*.csv;*.mat','Objective files (*.csv,*.mat)'}, ...
                     '目的行列(N×M)ファイルを選択');
if isequal(fn,0)
    fprintf('キャンセルしました。\n'); 
    return;
end

filepath = fullfile(fp, fn);
[~, base, ext] = fileparts(filepath);

% --- 目的行列を読み込み ---
switch lower(ext)
    case '.csv'
        F = readmatrix(filepath);
    case '.mat'
        S = load(filepath);
        cand = {'F','objs','ParetoObjs','pf_objs','gen_objs'};
        F = [];
        for k = 1:numel(cand)
            if isfield(S, cand{k}), F = S.(cand{k}); break; end
        end
        if isempty(F)
            error('MAT内に N×M の目的配列（F/objs/ParetoObjs/pf_objs/gen_objs）が見つかりません。');
        end
    otherwise
        error('未対応の拡張子: %s', ext);
end

% --- 目的名（列名） ---
obj_names = ["J_tau","J_pos","J_ori","J_sc","J_smooth"];

% --- セレクタ（評価基準） ---
% 既定：バランス3種＋単独最良（位置/姿勢）。好みに応じて default_selectors.m を編集。
selectors = default_selectors();

% --- 代表抽出 ---
[sel, T] = pick_representatives(F, ...
    'selectors', selectors, ...
    'obj_names', obj_names, ...
    'feasible_only', true, ...            % 可行(J_sc==0)のみで選抜
    'tol', 1e-12, ...
    'norm', 'percentile', 'plow',5, 'phigh',95);  % 外れ値に強い正規化

% --- 表示＆保存 ---
disp(T);
outcsv = fullfile(fp, [base '_representatives.csv']);
writetable(T, outcsv);
save(fullfile(fp, [base '_representatives.mat']), 'sel','T','obj_names','selectors');

fprintf('[OK] 代表一覧を保存: %s\n', outcsv);
assignin('base','T',T);
assignin('base','sel',sel);