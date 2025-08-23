function hv = hypervolumeContribution(F)

    [N, M] = size(F);
    
    ref_point = max(F, [], 1) + 1.0;

    hv_total = hypervolume(F, ref_point);

    hv = zeros(N, 1);
    for i = 1:N
        Fi = F;
        Fi(i, :) = [];
        hv_i = hypervolume(Fi, ref_point);
        hv(i) = hv_total - hv_i;
    end
end
