% 1) FK：壊れた入力は NaN で返る（クラッシュしない）
L = rand(1,6); th = [0 0 NaN 0 0 0];
[pos, q] = FK(L, th);  %#ok<NASGU>
assert(any(isnan(pos)));

% 2) quat_normalize：ゼロ四元数でも落ちない
qz = quaternion(0,0,0,0);
qn = quat_normalize(qz); %#ok<NASGU>

% 3) crowding_distance：NaN混在でも落ちず、無効は -Inf になる
F = [1 2; NaN 3; 4 5];
cd = crowding_distance(F);
assert(isinf(cd(2)) && cd(2)<0);

% 4) Evaluate：非有限は invalid_row に置換
% (あなたの Evaluate の呼び出し形に合わせて1行だけ試す)
if exist('Evaluate','file')
    pop = rand(2,6); pop(2,:) = NaN;
    J = Evaluate(pop);
    assert(all(isfinite(J(1,:))));
end