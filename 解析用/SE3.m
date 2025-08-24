classdef SE3
    methods (Static)
        function T = transl(x,y,z)
            T = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1];
        end

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

        function d = dotparts(q1, q2)
            a = compact(q1);
            b = compact(q2);
            d = dot(a, b);
        end
    end
end