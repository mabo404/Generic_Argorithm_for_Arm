classdef Kin
    methods (Static)
        function [p,q] = FK_safe(L, th)
            [p,q] = Kin.FK(L, th);
            p = real(double(p(:)));
            q = Quat.quat_normalize(q);
        end

        function [p,q] = FK_posori(L, th)
            [p,q] = Kin.FK_safe(L, th);
        end

        function J = numJacobian(L, th, f, q_des, Lc)
            epsFD = 1e-5;
            Jp = zeros(3,6); Jo = zeros(3,6);
            q_des = Quat.quat_normalize(q_des);

            for j = 1:6
                delta = max(epsFD, sqrt(eps)*max(1,abs(th(j))));
                thp = th; thm = th;
                thp(j) = th(j) + delta;
                thm(j) = th(j) - delta;

                [pp, qp] = f(L, thp); qp = Quat.quat_normalize(qp);
                [pm, qm] = f(L, thm); qm = Quat.quat_normalize(qm);

                Jp(:,j) = (pp - pm) / (2*delta) / Lc;
                wo_p = Quat.quat_log_err(q_des, qp);
                wo_m = Quat.quat_log_err(q_des, qm);
                Jo(:,j) = (wo_p - wo_m) / (2*delta);
            end
            J = [Jp; Jo];
        end

        function [pos, q_end] = FK(L, th)
            th = wrapToPi(real(double(th)));

            T = eye(4);
            for i = 1:6
                if ismember(i,[1,4,6])
                    R = SE3.rotz(th(i));
                else
                    R = SE3.rotx(th(i));
                end
                R = real(double(R));

                T = T * R * [eye(3) [0;0;L(i)]; 0 0 0 1];
            end

            pos   = T(1:3,4);
            R_end = real(double(T(1:3,1:3)));
            [U,~,V] = svd(R_end); R_end = U*V';
            if det(R_end) < 0, U(:,3) = -U(:,3); R_end = U*V'; end

            q_end = Quat.rotm2quat_safe(R_end);
        end

        function joints = FK_with_joints(L, theta)
            theta = real(double(theta(:)));
            if numel(L)~=6 || numel(theta)~=6
                joints = nan(3,7); return;
            end

            T = eye(4);
            joints = zeros(3,7);
            joints(:,1) = T(1:3,4);
            for i = 1:6
                if ismember(i,[1,4,6]), R = SE3.rotz(theta(i)); else, R = SE3.rotx(theta(i)); end
                T = T * R * SE3.transl(0,0,L(i));
                joints(:,i+1) = T(1:3,4);
            end

        end

        function [T_end, pos] = FK_with_pose(L, theta)
            L = real(double(L(:)));    theta = real(double(theta(:)));
            if numel(L)~=6 || numel(theta)~=6
                T_end = nan(4); pos = [NaN;NaN;NaN]; return;
            end

            T = eye(4);
            for i = 1:6
                if ismember(i, [1,4,6])
                    R = SE3.rotz(theta(i));
                else
                    R = SE3.rotx(theta(i));
                end
                T = T * R * SE3.transl(0, 0, L(i));
            end
            T_end = T;
            pos = T(1:3,4);
        end
    end
end