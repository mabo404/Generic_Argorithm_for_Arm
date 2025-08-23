clear; clc; close all;

target_gens = 28;
target_idx = 48;
logdir = 'run_20250815_165244';
rank_method = 'utopia';
animate = true;
L6_const = 0.1;

dt = 0.02; T_total = 10; t = 0:dt:T_total; r = t / T_total;
S = 10*(r.^3) - 15*(r.^4) + 6*(r.^5);
M = numel(S);

Xe = [-0.5; -0.5; 0.6];
QE = quaternion([-45 0 90], 'eulerd','ZYX','frame');
theta0 = deg2rad([-45 -90 90 0 90 0]);

for g = target_gens
    [x_best, objs_best] = load_individual_by_index(logdir, g, target_idx);
    L = [x_best(1:5).', L6_const];

    if ~(util_isfinite(L) && util_isfinite(theta0))
        warning('L/theta0 に非有限が含まれます。フレーム処理をスキップします。');
        continue;
    end

    x0 = [-(L(5) + L(6) - L(2)) * cos(pi/4);
          -(L(5) + L(6) - L(2)) * cos(pi/4);
           L(1) + L(3) + L(4)];
    Qs = quaternion([-45 0 90], 'eulerd','ZYX','frame');
    Qtraj = slerp(Qs, QE, S);
    Qtraj = arrayfun(@quat_normalize, Qtraj);
    Qtraj = Qtraj(:);
    qc = compact(Qtraj);
    for k = 2:M
        if dot(qc(k,:), qc(k-1,:)) < 0
            Qtraj(k) = quaternion(-qc(k,:));
            qc(k,:)  = -qc(k,:);
        end
    end

    qc = compact(Qtraj);
    bad = any(~isfinite(qc), 2);
    if any(bad)
        warning('Qtrajに非有限が含まれたため、該当点の姿勢を Qs に置換します。');
        Qtraj(bad) = Qs;
    end

    Pmid = 0.5*(x0 + Xe);
    Pmid(3) = Pmid(3) + 0.4;
    Ptraj = (1 - S).^2 .* x0 + 2*(1 - S).*S .* Pmid + S.^2 .* Xe;
    Ptraj(:,1) = x0;
    Qtraj(1) = Qs;
    
    if ~util_isfinite(Ptraj,'mat')
        error('Ptraj に NaN/Inf が含まれています。');
    end
    if ~all(arrayfun(@(q) util_isfinite(q,'quat'), Qtraj))
        error('Qtraj に NaN/Inf の四元数が含まれています。');
    end


    [theta, ~] = IK_solver(L, Ptraj, Qtraj, theta0);
    theta = unwrap(theta, [], 2);

    badcol = any(~isfinite(theta), 1);
    if any(badcol)
        warnN = nnz(badcol);
        warning('%d フレームで theta に NaN/Inf。該当フレームを可視化から除外します。', warnN);
        keep = ~badcol;
        Ptraj = Ptraj(:, keep);
        Qtraj = Qtraj(keep);
        theta = theta(:, keep);
    end
    M = size(theta, 2);
    if M == 0
        warning('有効なフレームがありません。処理を終了します。');
        return;
    end

    fig = figure('Name', sprintf('Generation %03d', g));
    set(fig, 'Position', [100 100 960 720]);
    reach = sum(L); margin = 0.5;
    obj_text = sprintf(' %.4g', objs_best); obj_text = ['[', strtrim(obj_text), ' ]'];
    annotation(fig, 'textbox', [0.02 0.83 0.45 0.14], ...
        'String', sprintf('Gen %d | obj = %s', g, obj_text), ...
        'FitBoxToText','on','BackgroundColor','w');

    if animate
    mp4name = sprintf('arm_gen_%03d_best.mp4', g);
    v = VideoWriter(mp4name, 'MPEG-4');
    v.FrameRate = 1/dt;
    open(v);
    end

    for i = 1:M

        if ~(util_isfinite(theta(:,i)) && util_isfinite(Ptraj(:,i),'vec3') && util_isfinite(Qtraj(i),'quat'))
            warning('frame %d をスキップ (non-finite input)', i);
            continue;
        end

        joints = FK_with_joints(L, theta(:,i));

        if ~util_isfinite(joints,'mat')
            warning('frame %d: joints が非有限のためスキップ', i);
            continue;
        end

        plot3(Ptraj(1,:), Ptraj(2,:), Ptraj(3,:), 'r--'); hold on;
        plot3(joints(1,:), joints(2,:), joints(3,:), 'bo-', 'LineWidth', 2);
        plot3(joints(1,end), joints(2,end), joints(3,end), 'go', 'MarkerSize', 8);

        [T_end, p_end] = FK_with_pose(L, theta(:,i));

        if ~(util_isfinite(T_end,'mat') && util_isfinite(p_end,'vec3'))
            warning('frame %d: FK_with_pose が非有限のためスキップ', i);
            continue;
        end

        R = T_end(1:3,1:3); sc = 0.1;
        quiver3(p_end(1), p_end(2), p_end(3), R(1,1), R(2,1), R(3,1), sc, 'r', 'LineWidth', 2);
        quiver3(p_end(1), p_end(2), p_end(3), R(1,2), R(2,2), R(3,2), sc, 'g', 'LineWidth', 2);
        quiver3(p_end(1), p_end(2), p_end(3), R(1,3), R(2,3), R(3,3), sc, 'b', 'LineWidth', 2);

        axis equal; grid on; view(37.5, 30);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        xlim([-reach-margin, reach+margin]);
        ylim([-reach-margin, reach+margin]);
        zlim([-reach-margin, reach+margin]);
        title(sprintf('Gen %d | step %d / %d', g, i, M));
        drawnow;

        if animate
            fr = getframe(fig);
            writeVideo(v, fr);
        else
            if i==M
                exportgraphics(gcf, sprintf('robot_gen_%03d.png', g), 'Resolution', 200);
            end
        end
        hold off;
    end

    if animate
        close(v);
    end
end

function [x_best, objs_best] = load_best_individual(logdir, g, rank_method)
    V = readmatrix(fullfile(logdir, sprintf('generation_%03d_vars.csv', g)));
    O = readmatrix(fullfile(logdir, sprintf('generation_%03d_objs.csv', g)));
    switch lower(rank_method)
        case 'obj1'
            [~, idx] = min(O(:,1));
        otherwise
            On = O;
            for j = 1:size(O,2)
                col = O(:,j); rng = max(col) - min(col);
                if rng < 1e-12, On(:,j) = 0; else, On(:,j) = (col - min(col))/rng; end
            end
            [~, idx] = min(vecnorm(On,2,2));
    end
    x_best   = V(idx,:).';
    objs_best= O(idx,:).';
end

function [x, objs] = load_individual_by_index(logdir, g, idx)
    V = readmatrix(fullfile(logdir, sprintf('generation_%03d_vars.csv', g)));
    O = readmatrix(fullfile(logdir, sprintf('generation_%03d_objs.csv', g)));
    assert(idx >= 1 && idx <= size(V,1), ...
        'target_idx=%d が範囲外です（1〜%d）。', idx, size(V,1));
    x    = V(idx,:).';
    objs = O(idx,:).';
end