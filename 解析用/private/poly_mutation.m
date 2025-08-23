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