clear; clc; close all;

target_gens = 49;
target_idx = 36;
logdir = 'run_20250826_113128';
animate = true;
L6_const = 0.1;

dt = 0.02;
T_total = 10;
T = 0:dt:T_total;
Xe = [-0.3; -0.3; 0.2];
r = T / T_total;
S = 10.*(r.^3) - 15.*(r.^4) + 6.*(r.^5);
Qs = quaternion([-45 0 90], 'eulerd','ZYX','frame');
QE = quaternion([-45 0 90], 'eulerd','ZYX','frame');

% 終点との符号合わせ（antipodal 対策）
dot_q = sum(compact(Qs).*compact(QE));
if dot_q < 0
    QE    = quaternion(-compact(QE));
    dot_q = -dot_q;
end

% slerp 特異域に近いときは nlerp、通常は slerp
if dot_q > 0.9995
    q1 = compact(Qs); q2 = compact(QE);      % 1x4 数値
    W1 = (1 - S(:)); W2 = S(:);              % Mx1
    Qnum = W1.*q1 + W2.*q2;                  % Mx4
    Qnum = Qnum ./ vecnorm(Qnum,2,2);
    Qtraj = quaternion(Qnum);                % 数値→quaternion配列
else
    Qtraj = slerp(Qs, QE, S);                % quaternionのまま補間
end
Qtraj = Qtraj(:);

% 連続性（符号）をそろえる
for k = 2:numel(Qtraj)
    if sum(compact(Qtraj(k)).*compact(Qtraj(k-1))) < 0
    Qtraj(k) = quaternion(-compact(Qtraj(k)));
    end
end

theta0 = deg2rad([-45 -90 90 0 90 0]);

for g = target_gens
    [x_best, objs_best] = load_individual_by_index(logdir, g, target_idx);
    L = [x_best(1:5).', L6_const];

    x0 = [-(L(5) + L(6) - L(2)) * cos(pi/4);
          -(L(5) + L(6) - L(2)) * cos(pi/4);
           L(1) + L(3) + L(4)];

    Pmid = 0.5*(x0 + Xe);
    Pmid(3) = Pmid(3) + 0.4;
    Ptraj = (1 - S).^2 .* x0 + 2*(1 - S).*S .* Pmid + S.^2 .* Xe;
    M = size(Ptraj, 2);
    [theta, ~] = IK.IK_solver(L, Ptraj, Qtraj(:), theta0);
    theta = unwrap(theta, [], 2);
    if isempty(theta) || size(theta,1)~=6 || size(theta,2)~=numel(T) || ~all(isfinite(theta(:)))
        error('IK結果が評価時の期待と不一致/非有限（サイズやNaN/Inf）。');
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
        joints = Kin.FK_with_joints(L, theta(:,i));

        plot3(Ptraj(1,:), Ptraj(2,:), Ptraj(3,:), 'r--'); hold on;
        plot3(joints(1,:), joints(2,:), joints(3,:), 'bo-', 'LineWidth', 2);
        plot3(joints(1,end), joints(2,end), joints(3,end), 'go', 'MarkerSize', 8);

        [T_end, p_end] = Kin.FK_with_pose(L, theta(:,i));

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
                exportgraphics(gcf, sprintf('robot_gen_best_%03d.png', g), 'Resolution', 200);
            end
        end
        hold off;
    end
    if animate
        close(v);
    end
end

function [x, objs] = load_individual_by_index(logdir, g, idx)
    V = readmatrix(fullfile(logdir, sprintf('gen_%04d_vars.csv', g)));
    O = readmatrix(fullfile(logdir, sprintf('gen_%04d_objs.csv', g)));
    assert(idx >= 1 && idx <= size(V,1), ...
        'target_idx=%d が範囲外です（1〜%d）。', idx, size(V,1));
    x    = V(idx,:).';
    objs = O(idx,:).';
end