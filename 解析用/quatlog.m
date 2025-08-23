function v = quatlog(q)
    if ~isa(q,'quaternion')
        q = quaternion(q(1),q(2),q(3),q(4));
    end
    q = quat_normalize(q);
    axang = quat2axang(q);
    axis  = axang(1:3)'; angle = axang(4);
    if angle < 1e-8, v = zeros(3,1); else, v = axis * angle; end
end
