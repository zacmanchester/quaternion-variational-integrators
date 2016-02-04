function q = qmult(q1, q2)
%QMULT performs quaternion multiplication. It follows the same order
%convention as matrix multiplication

v1 = q1(1:3);
s1 = q1(4);

v2 = q2(1:3);
s2 = q2(4);

q = [s1*v2 + s2*v1 + hat(v1)*v2; s1*s2 - v1'*v2];

end

