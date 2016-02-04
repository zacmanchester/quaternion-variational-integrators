function xdot = gyrostatODE(I, wheelfun, torquefun, t, x)

q = x(1:4);
qn = norm(q);
q = q/qn;
w = x(5:7);
r = x(8:10); %rotor momentum

%Get control torque
tau_r = wheelfun(t,x);
tau_ext = torquefun(t,x);

%Calculate omega derivative from Euler's equation + torques
wdot = -I\(hat(w)*(I*w+r) + tau_r - tau_ext);

%Calculate quaternion derivative from omega
qdot = .5*[hat(q(1:3)) + q(4)*eye(3); -q(1:3)']*w; %- (qn-1)*q;

xdot = [qdot; wdot; tau_r];

end

