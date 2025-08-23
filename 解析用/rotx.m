function R = rotx(th)
    R = [1 0       0        0;
         0 cos(th) -sin(th) 0;
         0 sin(th) cos(th)  0;
         0 0       0        1];
end