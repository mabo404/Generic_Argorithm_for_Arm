function q = quat_normalize(q)
    if isa(q,'quaternion')
        [a,b,c,d] = parts(q);
        n = sqrt(a.^2 + b.^2 + c.^2 + d.^2);
        n(n==0) = 1;
        a = a./n; b = b./n; c = c./n; d = d./n;
        q = quaternion(a,b,c,d);
    else
        v = double(q);
        if isvector(v) && numel(v)==4
            n = norm(v);
            if n == 0
                q = quaternion(1,0,0,0);
                return;
            end
            v = v./n;
            q = quaternion(v(1),v(2),v(3),v(4));
        else
            error('Unsupported quaternion format');
        end
    end
end
