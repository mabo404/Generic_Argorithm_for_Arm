function visualize_arm_pose(L, theta)
    if numel(L) ~= 6 || numel(theta) ~= 6
        error('L and theta must both be 6 elements');
    end

    joints = Kin.FK_with_joints(L, theta);

    T_end = eye(4);
    for i = 1:6
        if ismember(i, [1,4,6])
            R = rotz(theta(i));
        else
            R = rotx(theta(i));
        end
        T_end = T_end * R * transl(0, 0, L(i));
    end
    pos_end = T_end(1:3, 4);
    R_end = T_end(1:3, 1:3);

    figure; view(3); hold on;
    plot3(joints(1,:), joints(2,:), joints(3,:), 'bo-', 'LineWidth', 2);
    plot3(joints(1,end), joints(2,end), joints(3,end), 'go', 'MarkerSize', 8);

    scale = 0.1;
    quiver3(pos_end(1), pos_end(2), pos_end(3), R_end(1,1), R_end(2,1), R_end(3,1), scale, 'r', 'LineWidth', 2);
    quiver3(pos_end(1), pos_end(2), pos_end(3), R_end(1,2), R_end(2,2), R_end(3,2), scale, 'g', 'LineWidth', 2);
    quiver3(pos_end(1), pos_end(2), pos_end(3), R_end(1,3), R_end(2,3), R_end(3,3), scale, 'b', 'LineWidth', 2);

    axis equal; grid on;
    reach = sum(L);
    margin = 0.5;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    xlim([-reach-margin, reach+margin]);
    ylim([-reach-margin, reach+margin]);
    zlim([0, reach+margin]);
    title('Arm Pose Visualization');

    fprintf('End-effector position: X = %.3f, Y = %.3f, Z = %.3f\n', ...
        pos_end(1), pos_end(2), pos_end(3));

    q_end = rotm2quat(R_end);  % [w x y z] 形式（MATLAB標準）
    fprintf('End-effector orientation (quaternion): [w=%.3f, x=%.3f, y=%.3f, z=%.3f]\n', ...
        q_end(1), q_end(2), q_end(3), q_end(4));
end

L = [0.3 0.4 0.3 0.2 0.2 0.1];
theta = deg2rad([-45 -90 90 0 90 0])';

visualize_arm_pose(L, theta);