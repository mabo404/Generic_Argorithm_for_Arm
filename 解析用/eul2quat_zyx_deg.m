function q = eul2quat_zyx_deg(eul_deg)
z=eul_deg(1)*pi/180; y=eul_deg(2)*pi/180; x=eul_deg(3)*pi/180;
cz=cos(z/2); sz=sin(z/2); cy=cos(y/2); sy=sin(y/2); cx=cos(x/2); sx=sin(x/2);
w = cz*cy*cx + sz*sy*sx;
xq= cz*cy*sx - sz*sy*cx;
yq= cz*sy*cx + sz*cy*sx;
zq= sz*cy*cx - cz*sy*sx;
q = quat_normalize([w xq yq zq]);
end
