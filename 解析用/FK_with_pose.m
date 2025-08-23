function [T_end, pos] = FK_with_pose(L, theta)
    if ~(util_isfinite(L) && util_isfinite(theta))
        T_end = nan(4); pos = [NaN;NaN;NaN]; return;
    end

    L = real(double(L(:)));    theta = real(double(theta(:)));
    if numel(L)~=6 || numel(theta)~=6
        T_end = nan(4); pos = [NaN;NaN;NaN]; return;
    end

    T = eye(4);
    for i = 1:6
        if ismember(i, [1,4,6])
            R = rotz(theta(i));
        else
            R = rotx(theta(i));
        end
        T = T * R * transl(0, 0, L(i));
        if ~util_isfinite(T,'mat') || ~isreal(T)
            T_end = nan(4); pos = [NaN;NaN;NaN]; return;
        end
    end
    T_end = T;
    pos = T(1:3,4);
end