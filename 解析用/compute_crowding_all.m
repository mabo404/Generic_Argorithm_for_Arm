function cd_all = compute_crowding_all(F, rank)
    cd_all = zeros(size(F,1),1);
    for fl = 1:max(rank)
        idx = find(rank==fl);
        if ~isempty(idx)
            cd_all(idx) = crowding_distance(F(idx,:));
        end
    end
end