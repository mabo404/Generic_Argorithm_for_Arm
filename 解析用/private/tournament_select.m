function sel = tournament_select(tmp, rank, cd_all)
    n = size(tmp,1);
    sel = zeros(n,1);
    for i = 1:n
        cands = tmp(i,:);
        r = rank(cands);
        best = cands(r==min(r));           % ランク最小を抽出
        if numel(best)==1
            sel(i) = best;
        else
            [~,ix] = max(cd_all(best));     % 同率なら距離最大
            sel(i) = best(ix);
        end
    end
end