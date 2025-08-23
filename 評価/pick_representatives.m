function [sel, T] = pick_representatives(F, varargin)
%PICK_REPRESENTATIVES 多目的（最小化）結果から複数基準で代表個体を抽出
%
% [sel, T] = pick_representatives(F, 'Name',Value,...)
%
% 入力:
%   F            : N×M 目的行列（最小化）。既定列名:
%                  ["J_tau","J_pos","J_ori","J_sc","J_smooth"]
%
% Name-Value:
%   'selectors'      : struct配列（複数基準; default_selectors() 参照）
%   'obj_names'      : 1×M string 配列（目的名）
%   'tol'            : 可行判定 |J_sc|<=tol（既定 1e-12）
%   'feasible_only'  : true→可行集合のみから選抜（既定 true）
%   'norm'           : 'percentile'(既定) | 'minmax' | 'spec'
%   'plow','phigh'   : percentileの下上限（既定 5,95）
%   'spec_range'     : norm='spec' 時の [M×2]（各目的 [low high]）
%
% 出力:
%   sel              : sel.by_label.<label> に各基準の行インデックス
%   T                : 選ばれたユニーク個体の目的値表（index, tag 付き）

% ---------- オプション ----------
p = inputParser;
addParameter(p,'selectors',struct([]));
addParameter(p,'obj_names',["J_tau","J_pos","J_ori","J_sc","J_smooth"]);
addParameter(p,'tol',1e-12);
addParameter(p,'feasible_only',true);
addParameter(p,'norm','percentile');   % 'percentile'|'minmax'|'spec'
addParameter(p,'plow',5); addParameter(p,'phigh',95);
addParameter(p,'spec_range',[]);
parse(p,varargin{:});
opts = p.Results;

% ---------- 入力整形 ----------
[N,M] = size(F);
obj_names = opts.obj_names;
if numel(obj_names) ~= M, obj_names = "F" + string(1:M); end

% ---------- 可行集合（J_sc==0） ----------
colJsc = find(strcmp(obj_names,"J_sc"), 1);
if isempty(colJsc)
    error('obj_names に "J_sc" が見つかりません。');
end
feas_mask = abs(F(:,colJsc)) <= opts.tol;
if opts.feasible_only && any(feas_mask)
    used_idx = find(feas_mask);
else
    used_idx = (1:N).';
end
U  = F(used_idx,:);    % 実値（評価対象）

% ---------- 正規化 ----------
switch lower(opts.norm)
    case 'percentile'
        lo = prctile(U, opts.plow,  1);
        hi = prctile(U, opts.phigh, 1);
    case 'minmax'
        lo = min(U,[],1); hi = max(U,[],1);
    case 'spec'
        if isempty(opts.spec_range) || size(opts.spec_range,2)~=2 || size(opts.spec_range,1)~=M
            error('spec_range は [M×2] の [low high] を指定してください。');
        end
        lo = opts.spec_range(:,1).'; hi = opts.spec_range(:,2).';
    otherwise
        error('norm は ''percentile''|''minmax''|''spec'' のいずれかです。');
end
den = max(hi - lo, eps);
Un  = (U - lo) ./ den;
Un  = max(0, min(1.5, Un));   % 0..1.5 にクリップ（外れ値耐性）

% ---------- セレクタ ----------
selectors = opts.selectors;
if isempty(selectors)
    selectors = default_selectors();
end

chosen = []; labels = strings(0,1);
sel.by_label = struct();

for s = 1:numel(selectors)
    S = selectors(s);
    if ~isfield(S,'label') || strlength(string(S.label))==0
        S.label = S.type;
    end
    lbl = string(S.label);

    switch lower(S.type)
        case 'knee'
            idx = used_idx( argmin_vecnorm(Un,2) );

        case 'minsum'
            [~,rel] = min(sum(Un,2));
            idx = used_idx(rel);

        case 'minimax'
            [~,rel] = min(max(Un,[],2));
            idx = used_idx(rel);

        case 'objective_min'
            must_have(S,'objective','objective_min');
            j = name2idx(S.objective, obj_names);
            [~,rel] = min(U(:,j));
            cand = used_idx(U(:,j) == U(rel,j));
            idx  = tie_break(cand, Un, used_idx, getdef(S,'tie_breaker','knee'));

        case 'weighted_sum'
            must_have(S,'weights','weighted_sum');
            w = S.weights(:).';
            if numel(w) ~= M, error('weights の長さが列数 M と一致しません。'); end
            score = Un * w(:);
            [~,rel] = min(score);
            cand = used_idx(score == score(rel));
            idx  = tie_break(cand, Un, used_idx, getdef(S,'tie_breaker','knee'));

        case 'p_norm'
            pval = getdef(S,'p',2);
            w    = getdef(S,'weights',ones(1,M));
            if numel(w) ~= M, error('p_norm: weights の長さが M と一致しません。'); end
            WW = Un .* w;
            score = sum(abs(WW).^pval,2).^(1/pval);
            [~,rel] = min(score);
            cand = used_idx(score == score(rel));
            idx  = tie_break(cand, Un, used_idx, getdef(S,'tie_breaker','knee'));

        case 'lexi'
            must_have(S,'order','lexi');
            ord = arrayfun(@(x) name2idx(x, obj_names), S.order);
            keep = 1:size(U,1);
            for jj = 1:numel(ord)
                j = ord(jj);
                v = U(keep,j); m = min(v);
                keep = keep(v==m);
                if numel(keep)==1, break; end
            end
            cand = used_idx(keep);
            idx  = tie_break(cand, Un, used_idx, getdef(S,'tie_breaker','knee'));

        case 'threshold_then_knee'
            must_have(S,'thresholds','threshold_then_knee');
            thr = S.thresholds;
            if isstruct(thr)
                mask = true(size(U,1),1);
                fns = fieldnames(thr);
                for k = 1:numel(fns)
                    j = name2idx(fns{k}, obj_names);
                    if getdef(S,'use_normalized',false)
                        mask = mask & (Un(:,j) <= thr.(fns{k}));
                    else
                        mask = mask & (U(:,j)  <= thr.(fns{k}));
                    end
                end
            else
                if numel(thr) ~= M, error('thresholds は 1×M の数値ベクトルです。'); end
                if getdef(S,'use_normalized',false), mask = all(Un <= thr, 2);
                else,                                  mask = all(U  <= thr, 2);
                end
            end
            K = used_idx(mask);
            if isempty(K)
                idx = used_idx( argmin_vecnorm(Un,2) );      % 閾下が無ければ全体から膝
            else
                Krel = ismember(used_idx, K);
                UnK  = Un(Krel,:);
                pick = argmin_vecnorm(UnK,2);
                idx  = K(pick);
            end

        otherwise
            error('未知の selector.type: %s', S.type);
    end

    labels(end+1,1) = lbl;
    chosen(end+1,1) = idx;
    sel.by_label.(matlab.lang.makeValidName(lbl)) = idx;
end

% ---------- ユニーク化 & テーブル ----------
[uniq_idx, ~] = unique(chosen, 'stable');
vals = F(uniq_idx,:);
T = array2table(vals, 'VariableNames', cellstr(obj_names));
T.index = uniq_idx;
T = movevars(T, 'index', 'Before', 1);

tags = strings(numel(uniq_idx),1);
for i = 1:numel(uniq_idx)
    who = labels(chosen==uniq_idx(i));
    tags(i) = strjoin(who, ',');
end
T.tag = tags;

% ---------- 互換フィールド ----------
sel.knee_idx    = get_or_nan(sel.by_label, 'knee');
sel.minsum_idx  = get_or_nan(sel.by_label, 'minsum');
sel.minimax_idx = get_or_nan(sel.by_label, 'minimax');

end
