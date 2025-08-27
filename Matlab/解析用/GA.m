classdef GA
    methods (Static)
        function sel = tournament_select(tmp, rank, cd_all)
            n = size(tmp,1);
            sel = zeros(n,1);
            for i = 1:n
                cands = tmp(i,:);
                r = rank(cands);
                best = cands(r==min(r));
                if numel(best)==1
                    sel(i) = best;
                else
                    [~,ix] = max(cd_all(best));
                    sel(i) = best(ix);
                end
            end
        end

        function cd_all = compute_crowding_all(F, rank)
            cd_all = zeros(size(F,1),1);
            for fl = 1:max(rank)
                idx = find(rank==fl);
                if ~isempty(idx)
                    cd_all(idx) = crowding_distance(F(idx,:));
                end
            end
        end

        function cd = crowding_distance(F)
            N = size(F,1);
            M = size(F,2);

            cd = -Inf(N,1);
            row_ok = true(size(F,1), 1);
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

        function [next_vari, next_fit] = environmental_selection(combined_vari, combined_fit, pop_size)
            rank = pareto_rank(combined_fit);
            next_vari = []; next_fit = []; remain = pop_size; fl = 1;
            while remain > 0
                idx = find(rank==fl);
                if numel(idx) <= remain
                    next_vari = [next_vari; combined_vari(idx,:)];
                    next_fit  = [next_fit;  combined_fit(idx,:)];
                    remain = remain - numel(idx); fl = fl + 1;
                else
                    F_front = combined_fit(idx,:);
                    cd_local = crowding_distance(F_front);
                    [~,I] = sort(cd_local, 'descend');
                    sel = idx(I(1:remain));
                    next_vari = [next_vari; combined_vari(sel,:)];
                    next_fit  = [next_fit;  combined_fit(sel,:)];
                    break;
                end
            end
        end

        % function hv = hypervolume(F, ref_point)
        %     [N, M] = size(F);
        % 
        %     if M ~= 2 && M ~= 3
        %         error('hypervolume関数は2次元または3次元のみ対応しています');
        %     end
        % 
        %     F = sortrows(F);
        % 
        %     if M == 2
        %         hv = 0;
        %         prev_f1 = ref_point(1);
        %         for i = N:-1:1
        %             width = prev_f1 - F(i,1);
        %             height = ref_point(2) - F(i,2);
        %             hv = hv + width * height;
        %             prev_f1 = F(i,1);
        %         end
        %     else
        %         hv = 0;
        %         for i = 1:N
        %             vol = prod(ref_point - F(i,:));
        %             hv = hv + vol;
        %         end
        %     end
        % end
        % 
        % function hv = hypervolumeContribution(F)
        % 
        %     [N, M] = size(F);
        % 
        %     ref_point = max(F, [], 1) + 1.0;
        % 
        %     hv_total = hypervolume(F, ref_point);
        % 
        %     hv = zeros(N, 1);
        %     for i = 1:N
        %         Fi = F;
        %         Fi(i, :) = [];
        %         hv_i = hypervolume(Fi, ref_point);
        %         hv(i) = hv_total - hv_i;
        %     end
        % end

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

        function mutants = poly_mutation(pop, lb, ub, dis_m, pro_m)
            [n,v] = size(pop);
            rv  = rand(n,v);
            mu2 = rand(n,v);
            lb_full = repmat(lb, n, 1); ub_full = repmat(ub, n, 1);

            delta = min(pop - lb_full, ub_full - pop) ./ (ub_full - lb_full + eps);
            detaq = zeros(n,v);

            pos   = rv<=pro_m;
            left  = pos & (mu2<=0.5);
            right = pos & ~ (mu2<=0.5);
            detaq(left)  = (2*mu2(left) + (1-2*mu2(left)).*(1-delta(left)).^(dis_m+1)).^(1/(dis_m+1)) - 1;
            detaq(right) = 1 - (2*(1-mu2(right)) + 2*(mu2(right)-0.5).*(1-delta(right)).^(dis_m+1)).^(1/(dis_m+1));

            mutants = pop + detaq .* (ub_full - lb_full);
            mutants = min(max(mutants, lb_full), ub_full);
        end

        function children = sbx_crossover(parents, lb, ub, dis_c)
            n = size(parents,1); v = size(parents,2);
            assert(mod(n,2)==0,'parents must be even.');
            p1 = parents(1:2:end,:); p2 = parents(2:2:end,:);
            lb_half = repmat(lb, n/2, 1);
            mu1 = rand(n/2, v);

            beta  = (1 + 2.*min(p1,p2) - lb_half) ./ max(abs(p2 - p1), 1e-6);
            alpha = 2 - beta.^(-dis_c-1);
            inv_alpha = 1 ./ alpha;
            term1 = alpha .* mu1;
            term2 = 1 ./ (2 - alpha .* mu1);
            power_exp = 1/(dis_c+1);

            mask1 = (mu1 <= inv_alpha);
            mask2 = ~mask1;
            tmp1 = term1 .^ power_exp;
            tmp2 = term2 .^ power_exp;
            betaq = tmp1 .* mask1 + tmp2 .* mask2;
            betaq = betaq .* ((-1).^randi([0,1],n/2,v));

            c1 = 0.5*((1+betaq).*p1 + (1-betaq).*p2);
            c2 = 0.5*((1-betaq).*p1 + (1+betaq).*p2);
            children = [c1; c2];

            lb_full = repmat(lb, n, 1); ub_full = repmat(ub, n, 1);
            children = min(max(children, lb_full), ub_full);
        end

    end
end