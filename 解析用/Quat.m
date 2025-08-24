classdef Quat
    methods (Static)
        function v = quatlog(q)
            if isa(q,'quaternion')
                q = compact(Quat.quat_normalize(q));
            else
                q = Quat.quat_normalize(q);
                if isvector(q) && numel(q)==4, q = reshape(q,1,4); end
                if size(q,1) ~= 1 || size(q,2) ~= 4
                    error('quatlog:ExpectedSingle','Expected a single quaternion [1x4] or [4x1].');
                end
            end

            s   = max(-1, min(1, q(1)));
            v3  = q(2:4).';
            ang = 2*acos(s);

            if ang < 1e-8
                v = 2*v3;
            else
                denom = max(1e-12, sin(ang/2));
                axis  = v3/denom;
                v     = ang * axis;
                nv = norm(v);
                if nv > pi, v = v * (pi/nv); end
            end
        end

        function R = quat2rotm_safe(q)
            if isa(q,'quaternion'), q = compact(q); end
            q = double(q);
            if isvector(q) && numel(q)==4, q = reshape(q,1,4); end
            if size(q,2)~=4 && size(q,1)==4, q = q.'; end
            if size(q,2)~=4, error('quat2rotm_safe:UnsupportedInput','Expected [*,4].'); end

            N = size(q,1);
            if N > 1
                R = zeros(3,3,N);
                for i = 1:N
                    R(:,:,i) = Quat.quat2rotm_safe(q(i,:));  % 単発ケースに帰着
                end
            return;
            end

            q = Quat.quat_normalize(q);

            w = q(1); x = q(2); y = q(3); z = q(4);
            xx = x*x; yy = y*y; zz = z*z; xy = x*y; xz = x*z; yz = y*z; wx = w*x; wy = w*y; wz = w*z;

            R = [1-2*(yy+zz),   2*(xy - wz), 2*(xz + wy);
                2*(xy + wz), 1-2*(xx+zz),  2*(yz - wx);
                2*(xz - wy),  2*(yz + wx), 1-2*(xx+yy)];

            [U,~,V] = svd(R);
            R = U*V.';
            if det(R) < 0
                U(:,3) = -U(:,3);
                R = U*V.';
            end

            R = real(R);
        end

        function qs = quat_slerp(q0, q1, t)
            q0 = Quat.quat_normalize(q0); if isa(q0,'quaternion'), q0 = compact(q0); end
            q1 = Quat.quat_normalize(q1); if isa(q1,'quaternion'), q1 = compact(q1); end
            if dot(q0,q1) < 0, q1 = -q1; end
            t  = double(t);
            d  = max(-1,min(1,dot(q0,q1)));
            ang= acos(d);
            if ang < 1e-8, qs = repmat(q0,numel(t),1); return; end
            s0 = sin((1-t).*ang)/sin(ang);  s1 = sin(t.*ang)/sin(ang);
            s0 = s0(:);                     s1 = s1(:);
            qs = Quat.quat_normalize( s0.*q0 + s1.*q1 );
            if isa(qs,'quaternion'), qs = compact(qs); end
        end

        function q = quat_normalize(q)

            wasObj = isa(q,'quaternion');

            if wasObj
                [a,b,c,d] = parts(q);
                n = sqrt(a.^2 + b.^2 + c.^2 + d.^2);
                n(n==0) = 1;
                q = quaternion(a./n, b./n, c./n, d./n);
            else
                q = double(q);

                if isvector(q) && numel(q)==4, q = reshape(q,1,4); end
                if size(q,2)~=4 && size(q,1)==4, q = q.'; end
                if size(q,2)~=4
                    error('quat_normalize:UnsupportedInput','Expected [*,4] numeric.');
                end
                n = sqrt(sum(q.^2,2));
                n(n==0) = 1;
                q = q ./ n;
            end
        end

        function q = quat_mul(q1, q2)
            if isa(q1,'quaternion'), q1 = compact(q1); end
            if isa(q2,'quaternion'), q2 = compact(q2); end
            q1 = double(q1); q2 = double(q2);

            if isvector(q1) && numel(q1)==4, q1 = reshape(q1,1,4); end
            if isvector(q2) && numel(q2)==4, q2 = reshape(q2,1,4); end
            if size(q1,2)==1 && size(q1,1)==4, q1 = q1.'; end
            if size(q2,2)==1 && size(q2,1)==4, q2 = q2.'; end

            n1 = size(q1,1); n2 = size(q2,1);
            if n1 ~= n2
                if n1 == 1
                    q1 = repmat(q1, n2, 1);
                elseif n2 == 1
                    q2 = repmat(q2, n1, 1);
                else
                    error('quat_mul:sizeMismatch','Row counts differ: %d vs %d', n1, n2);
                end
            end

            w1 = q1(:,1); v1 = q1(:,2:4);
            w2 = q2(:,1); v2 = q2(:,2:4);
            w  = w1.*w2 - sum(v1.*v2, 2);
            v  = v1.*w2 + v2.*w1 + cross(v1, v2, 2);

            q = Quat.quat_normalize([w v]); 
        end

        function w = quat_log_err(q_des, q_cur)
            if isa(q_des,'quaternion'), q_des = compact(q_des); end
            if isa(q_cur,'quaternion'), q_cur = compact(q_cur); end
            q_des = Quat.quat_normalize(q_des);   % ← ここで [1x4] 数値に正規化
            q_cur = Quat.quat_normalize(q_cur);

            q_err = Quat.quat_mul(q_des, Quat.quat_conj(q_cur));
            q_err = Quat.quat_normalize(q_err);
            if isa(q_err,'quaternion'), q_err = compact(q_err); end

            v = q_err(2:4);
            s = max(-1, min(1, q_err(1)));
            ang = 2*acos(s);

            if ang < 1e-8
                w = 2*v(:);
            else
                axis = v / max(1e-12, sin(ang/2));
                w = ang * axis(:);
                if norm(w) > pi
                    w = w * (pi / norm(w));
                end
            end
        end

        function qc = quat_conj(q)
            if isa(q,'quaternion')
                q = compact(q);
            else
                q = double(q);
                if isvector(q) && numel(q)==4
                    q = reshape(q, 1, 4);
                elseif size(q,2)~=4 && size(q,1)==4
                    q = q.';
                end
            end
            qc = [q(:,1), -q(:,2:4)];
        end

        function q = rotm2quat_safe(R)
            R = real(double(R));

            [U,~,V] = svd(R);
            R = U*V.';
            if det(R) < 0
                U(:,3) = -U(:,3);
                R = U*V.';
            end
 
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
            q = Quat.quat_normalize([w x y z]);
        end

        function q = eul2quat_zyx_deg(eul_deg)
            e = double(eul_deg);
            if isvector(e) && numel(e)==3
                e = reshape(e, 1, 3);
            end
            if size(e,2) ~= 3
                error('eul2quat_zyx_deg:ExpectedNx3','Expected [*,3] Euler angles [Z Y X] in degrees.');
            end

            h = e * (pi/180/2);
            z = h(:,1); y = h(:,2); x = h(:,3);
            cz = cos(z);  sz = sin(z);
            cy = cos(y);  sy = sin(y);
            cx = cos(x);  sx = sin(x);

            w  = cz.*cy.*cx + sz.*sy.*sx;
            xq = cz.*cy.*sx - sz.*sy.*cx;
            yq = cz.*sy.*cx + sz.*cy.*sx;
            zq = sz.*cy.*cx - cz.*sy.*sx;

            q = Quat.quat_normalize([w xq yq zq]);
            if isa(q,'quaternion'), q = compact(q); end
        end
    end
end