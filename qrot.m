function rrot = qrot(q, r)
%QROT Rotates a 3D vector r by a quaternion q

v = q(1:3);
vh = hat(v);
w = q(4);

rrot = r + 2*vh*(vh*r + w*r);

end

