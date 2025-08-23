function [theta, info] = IK_solver(L, Ptraj, Qtraj, theta0)
    M = size(Ptraj,2);
    if size(Qtraj,2) ~= 1 || ~isa(Qtraj, 'quaternion')
        if isa(Qtraj,'quaternion') && size(Qtraj,1) == M
        elseif isa(Qtraj,'quaternion') && size(Qtraj,2) == M
            Qtraj = Qtraj.';
        else
            error('Qtraj must be Mx1 quaternion array');
        end
    end
    if ~(util_isfinite(L) && util_isfinite(Ptraj,'mat') && util_isfinite(theta0) && isa(Qtraj,'quaternion') && size(Qtraj,1)==M)
        info  = struct('status','badInput');
        theta = nan(6,M);
        return;
    end
    Qtraj = quat_normalize(Qtraj);

    theta = zeros(6,M);
    th    = real(double(theta0(:)));

    Lc = sum(L);
    maxIter = 120;
    tol_pos = 1e-4;
    tol_ori = 1e-5;
    maxStep = 30*pi/180;
    back_beta = 0.7;
    back_max = 12;
    lambda_min = 1e-5;
    k_lambda = 0.2;

    info = struct();
    info.iter = zeros(1,M);
    info.status_per_t = repmat({''},1,M);

    for t = 1:M
        p_des = Ptraj(:,t);
        q_des = Qtraj(t,1);

        th_best = th; F_best = inf; status = 'running';

        for it = 1:maxIter
            [p_cur, q_cur] = FK_safe(L, th);
            if ~(util_isfinite(p_cur,'vec3') && util_isfinite(q_cur,'quat'))
                status='fkNaN'; break;
            end


            e_p = (p_des - p_cur) / Lc;
            e_o = quatLogErr(q_des, q_cur);
            e = [e_p; e_o];

            if norm(e(1:3)) < tol_pos && norm(e(4:6)) < tol_ori
                status='ok';
                break;
            end

            J = numJacobian(L, th, @FK_posori, q_des, Lc);
            if ~util_isfinite(J,'mat')
                status='JNaN'; break;
            end

            w_pos = 1.0;
            w_ori = 0.05;

            W  = diag([w_pos w_pos w_pos w_ori w_ori w_ori]);
            eW = W * e;
            JW = W * J;

            lambda = max(lambda_min, k_lambda*norm(eW));
            A      = JW*JW' + (lambda^2)*eye(6);
            if rcond(A) < 1e-10
                lambda = max(lambda, 1e-1);
                A = J*J' + (lambda^2)*eye(6);
            end
            dth    = JW' * (A \ eW);

            scale = max(1, max(abs(dth)/maxStep));
            dth_limited = dth / scale;

            F0 = 0.5*norm(eW)^2; g = JW.'*eW; step = 1.0; accepted = false;
            th_cand = th; e_cand = e; F_cand = F0;
            for k = 1:back_max
                th_try = th + step*dth_limited;
                [p_try, q_try] = FK_posori(L, th_try);
                if ~(util_isfinite(p_try,'vec3') && util_isfinite(q_try,'quat'))
                    step = step * back_beta; continue;
                end

                e_try = [ (p_des - p_try)/Lc ; quatLogErr(q_des, q_try) ];
                eW_try = W * e_try;
                F1 = 0.5*norm(eW_try)^2;
                if F1 <= F0 - 1e-3*step*(norm(g)^2)
                    th_cand = th_try; e_cand = e_try; F_cand = F1; accepted = true; break;
                end
                step = step * back_beta;
            end

            if accepted
                th = th_cand; e = e_cand; F0 = F_cand;
                if F_cand < F_best && isfinite(F_cand)
                    th_best = th; F_best = F_cand;
                end
                lambda_min = max(1e-5, 0.7*lambda_min);
            else
                lambda_min = min(1.0, 2*lambda_min);
                if norm(g,inf) < 1e-6
                    status='no_progress'; break;
                end
            end

            if norm(J.'*e,inf) < 1e-5
                status='stationary'; break;
            end
            if norm(dth_limited,inf) < 1e-6
                status='stalled'; break;
            end
        end

        if isfinite(F_best)
            theta(:,t) = th_best;
        else
            theta(:,t) = th;
        end
        info.iter(t)        = it;
        info.status_per_t{t}= status;

        th = theta(:,t);

        if any(strcmp(status, {'fkNaN','JNaN'}))
            if t < M, theta(:,t+1:end) = repmat(theta(:,t), 1, M - t); end
            info.status = status; info.t = t; return;
        end
    end

    if ~isfield(info,'status'), info.status = 'ok'; end
    info.maxIter = maxIter;
end

function [p,q] = FK_safe(L, th)
    [p,q] = FK(L, th);
    p = real(double(p(:)));
    q = quat_normalize(q);
end

function [p,q] = FK_posori(L, th)
    [p,q] = FK_safe(L, th);
end

function J = numJacobian(L, th, f, q_des, Lc)
    epsFD = 1e-5;
    Jp = zeros(3,6); Jo = zeros(3,6);
    q_des = quat_normalize(q_des);

    for j = 1:6
        delta = max(epsFD, sqrt(eps)*max(1,abs(th(j))));
        thp = th; thm = th;
        thp(j) = th(j) + delta;
        thm(j) = th(j) - delta;

        [pp, qp] = f(L, thp); qp = quat_normalize(qp);
        [pm, qm] = f(L, thm); qm = quat_normalize(qm);

        Jp(:,j) = (pp - pm) / (2*delta) / Lc;
        wo_p = quatLogErr(q_des, qp);
        wo_m = quatLogErr(q_des, qm);
        Jo(:,j) = (wo_p - wo_m) / (2*delta);
    end
    J = [Jp; Jo];
end

function q = quat_normalize(q)
    if isa(q,'quaternion')
        [a,b,c,d] = parts(q);
        n = sqrt(a.^2 + b.^2 + c.^2 + d.^2);
        a = a./n; b = b./n; c = c./n; d = d./n;
        q = quaternion(a,b,c,d);
    else
        v = double(q);
        if isvector(v) && numel(v)==4
            n = norm(v);
            v = v./n; q = quaternion(v(1),v(2),v(3),v(4));
        else
            error('Unsupported quaternion format');
        end
    end
end

function w = quatLogErr(qd, qc)
    qd = quat_normalize(qd); qc = quat_normalize(qc);
    if qdot(qd,qc) < 0
    qc = quaternion(-compact(qc));
    end
    q = qd * conj(qc);
    [a,b,c,d] = parts(q);
    a = max(-1,min(1,a));
    ang = 2*acos(a);
    s = sqrt(max(1 - a.*a, 0));
    if s < 1e-12
        w = [0;0;0];
    else
        w = ang * [b;c;d] / s;
    end
end

function c = qdot(q1,q2)
    [a1,b1,c1,d1] = parts(q1);
    [a2,b2,c2,d2] = parts(q2);
    c = a1.*a2 + b1.*b2 + c1.*c2 + d1.*d2;
end