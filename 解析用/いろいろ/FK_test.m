clearvars; clc; close all;
function X = FK(L,theta)
    T = eye(4);
    for i=1:6
        if ismember(i, [1,4,6])
            R = rotz(theta(i));
        else
            R = rotx(theta(i));
        end
        T = T * R * transl(0,0,L(i));
    end
    X = T(1:3,4);
end
L = [0.3 0.3 0.3 0.3 0.3 0.1];
thetai = [0 pi/4 0 0 0 0]';
disp(FK(L, thetai));
theta = [0 -pi/4 0 0 0 0]';
disp(FK(L, theta));
function R = rotz(th)
    R = [cos(th) -sin(th) 0 0;
        sin(th)  cos(th)  0 0;
        0        0        1 0;
        0        0        0 1];
end
function R = rotx(th)
    R = [1 0       0        0;
         0 cos(th) -sin(th) 0;
         0 sin(th) cos(th)  0;
         0 0       0        1];
end
function T = transl(x,y,z)
    T = [1 0 0 x;
         0 1 0 y;
         0 0 1 z;
         0 0 0 1];
end