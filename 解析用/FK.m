function [pos, q_end] = FK(L, th)
    if ~(util_isfinite(L) && util_isfinite(th))
        pos = [NaN;NaN;NaN]; q_end = [NaN NaN NaN NaN]; return;
    end
    th = wrapToPi(real(double(th)));

    T = eye(4);
    for i = 1:6
        if ~util_isfinite(th) || any(abs(th) > 100*pi)
            pos = [NaN;NaN;NaN]; q_end = [NaN NaN NaN NaN]; return;
        end


        if ismember(i,[1,4,6])
            R = rotz(th(i));
        else
            R = rotx(th(i));
        end
        R = real(double(R));

        T = T * R * [eye(3) [0;0;L(i)]; 0 0 0 1];

        if ~util_isfinite(T,'mat') || ~isreal(T)
            pos = [NaN;NaN;NaN]; q_end = [NaN NaN NaN NaN]; return;
        end
    end

    pos   = T(1:3,4);
    R_end = real(double(T(1:3,1:3)));
    [U,~,V] = svd(R_end); R_end = U*V';
    if det(R_end) < 0, U(:,3) = -U(:,3); R_end = U*V'; end

    q_end = rotm2quat_safe(R_end);
end
