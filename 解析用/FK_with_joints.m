function joints = FK_with_joints(L, theta)
    if ~(util_isfinite(L) && util_isfinite(theta))
        joints = nan(3,7); return;
    end
    theta = real(double(theta(:)));
    if numel(L)~=6 || numel(theta)~=6
        joints = nan(3,7); return;
    end

    T = eye(4);
    joints = zeros(3,7);
    joints(:,1) = T(1:3,4);
    for i = 1:6
        if ismember(i,[1,4,6]), R = rotz(theta(i)); else, R = rotx(theta(i)); end
        T = T * R * transl(0,0,L(i));
        if ~util_isfinite(T,'mat') || ~isreal(T)
            joints(:,i+1:end) = NaN;   % 残りはNaNで埋めて脱出
            return;
        end
        joints(:,i+1) = T(1:3,4);
    end
    if ~util_isfinite(joints,'mat'), joints(:) = NaN; end
end
