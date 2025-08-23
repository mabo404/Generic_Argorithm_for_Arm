function qc = quat_conj(q)
if isa(q,'quaternion'), q = compact(q); end
qc = [q(:,1), -q(:,2:4)];
end