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
