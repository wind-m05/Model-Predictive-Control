function [outputArg1,outputArg2] = solve_ss(inputArg1,inputArg2)
%SOLVESS Summary of this function goes here
%   Detailed explanation goes here
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
for k = 1:k_sim
    uk(:,k) = K_mpc*xk(:,k);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

end

