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