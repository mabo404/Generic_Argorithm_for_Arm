function q = rotm2quat_safe(R)
    R = real(double(R));
    t = R(1,1)+R(2,2)+R(3,3);
    if t > 0
        s = 0.5 / sqrt(max(1e-12, t+1)); w = 0.25/s;
        x = (R(3,2)-R(2,3))*s; y = (R(1,3)-R(3,1))*s; z = (R(2,1)-R(1,2))*s;
    else
        [~,i] = max([R(1,1) R(2,2) R(3,3)]);
        switch i
        case 1
            s=2*sqrt(max(1e-12,1+R(1,1)-R(2,2)-R(3,3)));
            w=(R(3,2)-R(2,3))/s; x=0.25*s; y=(R(1,2)+R(2,1))/s; z=(R(1,3)+R(3,1))/s;
        case 2
            s=2*sqrt(max(1e-12,1-R(1,1)+R(2,2)-R(3,3)));
            w=(R(1,3)-R(3,1))/s; x=(R(1,2)+R(2,1))/s; y=0.25*s; z=(R(2,3)+R(3,2))/s;
        otherwise
            s=2*sqrt(max(1e-12,1-R(1,1)-R(2,2)+R(3,3)));
            w=(R(2,1)-R(1,2))/s; x=(R(1,3)+R(3,1))/s; y=(R(2,3)+R(3,2))/s; z=0.25*s;
        end
    end
    q = quat_normalize([w x y z]);
    if isa(q,'quaternion'), q = compact(q); end
end