function [t, qhist, whist, phihist] = GyrostatQVI(I, q0, w0, rhohist, tauhist, dt, tspan)
%GyrostatQVI is a quaternion variational integrator for the motion of a
%gyrostat (a rigid body with internal rotors/reaction wheels).
%   I = 3x3 Inertia matrix in body frame
%   q0 = [q1 q2 q3 q0]' initial attitude quaternion (body to inertial)
%   w0 = initial body angular velocity in body frame
%   rhohist = 3xlength(t) vector of internal rotor momentum in body frame
%   tauhist = 3xlength(t) vector of external torque inputs in the body frame
%   dt = timestep
%   tspan = [t0 tfinal] timespan of simulation

%Initialize Variables
t = tspan(1):dt:tspan(2);
qhist = zeros(4, length(t));
whist = zeros(3, length(t));
phihist = zeros(3, length(t)-1);
qhist(:,1) = q0;
whist(:,1) = w0;

Iinv = inv(I);
p = I*w0 + rhohist(:,1);

phi = w0*dt/2; %initial guess
for k = 1:(length(t)-1)
    
    %Use Newton's method to calculate phi
    g = .5*dt*p + .5*dt*dt*tauhist(:,k+1);
    for j = 1:3
        e = sqrt(1-phi'*phi)*(I*phi + (dt/2)*rhohist(:,k+1)) + hat(phi)*(I*phi + (dt/2)*rhohist(:,k+1)) - g;
        dedphi = sqrt(1-phi'*phi)*I - I*(phi*phi')/sqrt(1-phi'*phi) + hat(phi)*I - hat(I*phi) - .5*dt*rhohist(:,k+1)*phi'/sqrt(1-phi'*phi) - .5*dt*hat(rhohist(:,k));
        phi = phi - dedphi\e;
    end
    
    f = [phi; sqrt(1-phi'*phi)];
    qhist(:,k+1) = qmult(qhist(:,k),f);
    
    p = (2/dt)*(sqrt(1-phi'*phi)*(I*phi + (dt/2)*rhohist(:,k+1)) - hat(phi)*(I*phi + (dt/2)*rhohist(:,k+1)));
    whist(:,k+1) = Iinv*(p-rhohist(:,k+1));
    phihist(:,k) = phi;
    
end

end