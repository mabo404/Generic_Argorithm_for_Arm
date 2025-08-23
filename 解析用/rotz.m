function R = rotz(th)
    R = [cos(th) -sin(th) 0 0;
        sin(th)  cos(th)  0 0;
        0        0        1 0;
        0        0        0 1];
end