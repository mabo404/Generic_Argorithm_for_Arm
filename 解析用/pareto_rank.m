function rank = pareto_rank(F)
    N = size(F,1);
    rank = inf(N,1);
    dominated_count = zeros(N,1);
    dominates_list = cell(N,1);
    front{1} = [];
    for p=1:N
        Sp = [];
        np = 0;
        for q=1:N
            if all(F(p,:)<=F(q,:)) && any(F(p,:)<F(q,:))
                Sp = [Sp, q];
            elseif all(F(q,:)<=F(p,:)) && any(F(q,:)<F(p,:))
                np = np + 1;
            end
        end
        dominates_list{p} = Sp;
        dominated_count(p) = np;
        if np == 0
            rank(p) = 1;
            front{1} = [front{1}, p];
        end
    end
    i = 1;
    while ~isempty(front{i})
        Q = [];
        for p = front{i}
            for q = dominates_list{p}
                dominated_count(q) = dominated_count(q) - 1;
                if dominated_count(q) == 0
                    rank(q) = i + 1;
                    Q = [Q, q];
                end
            end
        end
        i = i + 1;
        front{i} = Q;
    end
end