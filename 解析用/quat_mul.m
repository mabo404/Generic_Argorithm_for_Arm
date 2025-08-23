function q = quat_mul(q1, q2)
    if isa(q1,'quaternion'), q1 = compact(q1); end
    if isa(q2,'quaternion'), q2 = compact(q2); end
    w1=q1(:,1); v1=q1(:,2:4);
    w2=q2(:,1); v2=q2(:,2:4);
    w  = w1.*w2 - sum(v1.*v2,2);
    v  = v1.*w2 + v2.*w1 + cross(v1,v2,2);
    q  = quat_normalize([w v]);
end