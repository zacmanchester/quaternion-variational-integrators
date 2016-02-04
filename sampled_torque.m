function tau = sampled_torque(tau_hist, dt, t, x)
%Implements a zero-order hold sampled torque input

sample = t/dt;

index = 1+floor(sample);

tau = tau_hist(:,index);

end

