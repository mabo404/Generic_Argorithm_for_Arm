function w = quat_log_err(q_des, q_cur)
    q_err = quat_mul(q_des, quat_conj(q_cur));
    q_err = quat_normalize(q_err);
    if isa(q_err,'quaternion'), q_err = compact(q_err); end
    v = q_err(2:4); s = max(-1,min(1,q_err(1)));
    ang = 2*acos(s);
    if ang < 1e-8
        w = [0;0;0];
    else
        axis = v / max(1e-12, sin(ang/2));
        w = ang * axis(:);
        if norm(w) > pi, w = w * (pi/norm(w)); end
    end
end