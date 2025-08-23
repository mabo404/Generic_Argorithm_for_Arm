function cd = crowding_distance(F)
    N = size(F,1);
    M = size(F,2);

    cd = -Inf(N,1);
    row_ok = all(isfinite(F), 2);
    if sum(row_ok) <= 2
        cd(row_ok) = Inf;
        return;
    end

    if N <= 2
        cd(:) = Inf;
        return;
    end

    for m = 1:M
        Iv = find(row_ok);
        [~, ord] = sort(F(Iv, m));
        I = Iv(ord);

        cd(I(1))  = Inf;
        cd(I(end)) = Inf;

        fmin = F(I(1), m);
        fmax = F(I(end), m);
        range = fmax - fmin;
        if range <= 0
            continue;
        end
        denom = range;

        for k = 2:numel(I)-1
            cd(I(k)) = cd(I(k)) + (F(I(k+1),m) - F(I(k-1),m)) / denom;
        end
    end
end
