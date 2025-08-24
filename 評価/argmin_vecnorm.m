function r = argmin_vecnorm(X,p)
%ARGMIN_VECNORM 行ベクトルの p-ノルム最小の行インデックス（1始まり）
    if nargin<2, p = 2; end
    nrm = sum(abs(X).^p,2).^(1/p);
    [~,r] = min(nrm);
end
