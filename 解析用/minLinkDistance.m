function dmin = minLinkDistance(L,theta)
    [~,N] = size(theta);
    dmin = inf;
    for k=1:N
        P = zeros(3,7);
        T4 = eye(4); 
        P(:,1)=[0;0;0];
        for i=1:6
            if ismember(i, [1,4,6])
                R = rotz(theta(i,k));
            else
                R = rotx(theta(i,k));
            end
            T4 = T4 * R * transl(0,0,L(i));
            P(:,i+1) = T4(1:3,4);
        end
        for i=1:6
            for j=i+2:6
                dmin = min(dmin, segmentSegmentDistance(P(:,i),P(:,i+1),P(:,j),P(:,j+1)));
            end
        end
    end
end