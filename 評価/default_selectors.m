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

seldefs = repmat(tpl, 1, 7);

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

seldefs(4).type       = 'objective_min';
seldefs(4).label      = 'best_smooth';
seldefs(4).objective  = 'J_smooth';
seldefs(4).tie_breaker= 'knee';

seldefs(5).type       = 'weighted_sum';
seldefs(5).label      = 'best_pos_ori';
seldefs(5).weights    = [0, 0.5, 0.5, 0, 0];

seldefs(6).type       = 'weighted_sum';
seldefs(6).label      = 'best_tau_pos_ori';
seldefs(6).weights    = [1/3, 1/3, 1/3, 0, 0];

seldefs(7).type       = 'p_norm';
seldefs(7).label      = 'best';
seldefs(7).p          = 2;
seldefs(7).weights    = [0.25, 0.25, 0.25, 0, 0.25];

end
