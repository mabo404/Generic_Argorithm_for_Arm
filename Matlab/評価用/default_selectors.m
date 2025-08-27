function seldefs = default_selectors()
tpl = struct( ...
    'type',           '', ...
    'label',          '', ...
    'objective',      [], ...
    'weights',        [], ...
    'p',              [], ...
    'order',          [], ...
    'thresholds',     [], ...
    'use_normalized', [], ...
    'tie_breaker',    ''  ...
);

seldefs = repmat(tpl, 1, 6);

seldefs(1).type       = 'objective_min';
seldefs(1).label      = 'best_tau';
seldefs(1).objective  = 'J_tau';
seldefs(1).tie_breaker= 'knee';

seldefs(2).type       = 'objective_min';
seldefs(2).label      = 'best_pos';
seldefs(2).objective  = 'J_pos';
seldefs(2).tie_breaker= 'knee';

seldefs(3).type       = 'objective_min';
seldefs(3).label      = 'best_ori';
seldefs(3).objective  = 'J_ori';
seldefs(3).tie_breaker= 'knee';

seldefs(4).type       = 'weighted_sum';
seldefs(4).label      = 'best_pos_ori';
seldefs(4).weights    = [0, 0.5, 0.5, 0];

seldefs(5).type       = 'weighted_sum';
seldefs(5).label      = 'best_tau_pos_ori';
seldefs(5).weights    = [1/3, 1/3, 1/3, 0];

seldefs(6).type       = 'p_norm';
seldefs(6).label      = 'best';
seldefs(6).p          = 2;
seldefs(6).weights    = [1/3, 1/3, 1/3, 0];

end