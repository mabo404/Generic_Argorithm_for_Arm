function idx = tie_break(candidates, Un, used_idx, how)
%TIE_BREAK 同値解の解消（'knee'|'minsum'|'minimax'）
    if numel(candidates) <= 1
        if isempty(candidates), idx = []; else, idx = candidates; end
        return;
    end
    cand_rel = ismember(used_idx, candidates);
    Cn = Un(cand_rel,:);
    switch lower(how)
        case 'minsum'
            [~,r] = min(sum(Cn,2));
        case 'minimax'
            [~,r] = min(max(Cn,[],2));
        otherwise % 'knee'
            r = argmin_vecnorm(Cn,2);
    end
    idx = candidates(r);
end
