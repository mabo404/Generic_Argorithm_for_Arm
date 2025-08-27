classdef Collision
    methods (Static) 
        function dmin = minLinkDistance(L,theta)
            [~,N] = size(theta);
            dmin = inf;
            for k=1:N
                P = zeros(3,7);
                T4 = eye(4); 
                P(:,1)=[0;0;0];
                for i=1:6
                    if ismember(i, [1,4,6])
                        R = SE3.rotz(theta(i,k));
                    else
                        R = SE3.rotx(theta(i,k));
                    end
                    T4 = T4 * R * SE3.transl(0,0,L(i));
                    P(:,i+1) = T4(1:3,4);
                end
                for i=1:6
                    for j=i+2:6
                        dmin = min(dmin, Collision.segmentSegmentDistance(P(:,i),P(:,i+1),P(:,j),P(:,j+1)));
                    end
                end
            end
        end

        function d = segmentSegmentDistance(A,B,C,D)
            u = B-A;
            v = D-C;
            w = A-C;
            a = dot(u,u);
            b = dot(u,v);
            c = dot(v,v);
            d_ = dot(u,w);
            e = dot(v,w);
            D0 = a*c - b^2;
            if D0<1e-8
                sc=0; tc=(b>c)*(d_/b)+(b<=c)*(e/c);
            else
                sc=(b*e - c*d_)/D0; tc=(a*e - b*d_)/D0;
            end
            sc = min(max(sc,0),1); tc = min(max(tc,0),1);
            dp = w + sc*u - tc*v;
            d = norm(dp);
        end
    end
end