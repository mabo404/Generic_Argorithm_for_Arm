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