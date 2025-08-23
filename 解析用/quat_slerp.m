function qs = quat_slerp(q0, q1, t)
    q0=quat_normalize(q0); if isa(q0,'quaternion'), q0 = compact(q0); end
    q1=quat_normalize(q1); if isa(q1,'quaternion'), q1 = compact(q1); end
    if dot(q0,q1) < 0, q1 = -q1; end
    d = max(-1,min(1,dot(q0,q1)));
    ang = acos(d);
    if ang < 1e-8, qs = repmat(q0,numel(t),1); return; end
    s0 = sin((1-t)*ang)/sin(ang); s1 = sin(t*ang)/sin(ang);
    qs = quat_normalize( s0.*q0 + s1.*q1 );
    if isa(qs,'quaternion'), qs = compact(qs); end
end