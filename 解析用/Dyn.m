classdef Dyn
    methods (Static)
        function torque = computeTorque(L, theta, ddtheta, m_joint)
            [n, M] = size(theta);
            torque = zeros(n, M);

            D   = 0.04;
            t   = 0.005;
            d   = D - 2*t;
            rho = 2700;
            A   = (pi/4) * (D^2 - d^2);

            L   = L(:);
            m_link = rho * A .* L;
            m_all  = m_link + m_joint(:);

            c = 0.5 * ones(1, n);

            Ivals = 0.5 * m_all .* ((D^2 + d^2) / 4);
            Ivals = max(Ivals, 1e-6);

            gvec  = [0; 0; -9.81];
            epsFD = 1e-5;

            for k = 1:M
                tau_g = zeros(n,1);

                for ell = 1:n
                    Jpos = zeros(3, ell);

                    for j = 1:ell
                        tp = theta(:,k); tp(j) = tp(j) + epsFD;
                        tm = theta(:,k); tm(j) = tm(j) - epsFD;

                        pp = Dyn.fk_link_cog(L, tp, ell, c);
                        pm = Dyn.fk_link_cog(L, tm, ell, c);

                        Jpos(:, j) = (pp - pm) / (2*epsFD);
                    end

                    F_ell = m_all(ell) * gvec;
                    tau_g(1:ell) = tau_g(1:ell) + Jpos.' * F_ell;
                end

                tau_iner     = Ivals .* ddtheta(:,k);
                torque(:,k)  = tau_iner + tau_g;
            end

            if any(isnan(torque(:)))
                warning('torqueにNaNが含まれています');
            end
            if any(~isreal(torque(:)))
                warning('torqueに虚数が含まれています');
            end
        end

        function p = fk_link_cog(L, theta, ell, c)
            T = eye(4);
            for i = 1:ell
                if ismember(i, [1,4,6])
                    R = rotz(theta(i));
                else
                    R = rotx(theta(i));
                end
                T = T * R * transl(0,0,L(i));
            end
            ph = T * [0; 0; c(ell)*L(ell); 1];
            p  = ph(1:3);
        end
    end
end