% 5LMB0 Model predictive control
% Instruction 1
clear all
close all

%% 3.1 Plant model
% (a) State-space representation: Select one model
% First-order model
    Ac = 0;
    Bc = 1;
    Cc = 1;
    Dc = 0;
    x0 = [1]; % Initial point
% Second-order model
    Ac = [0 1;0 0];
    Bc = [0;-1];
    Cc = eye(2);
    Dc = [0;0];
    x0 = [3;1]; % Initial point

% (b) Discretization
Ts = 0.1;
ss_c = ss(Ac,Bc,Cc,Dc);
ss_d = c2d(ss_c,Ts,'zoh');
A = ss_d.A;
B = ss_d.B;
C = ss_d.C;
D = ss_d.D;
[nx,nu] = size(B);

% (c) Controllability
ctrb_M = ctrb(A,B); % check if ctrb_M is full rank

% (d) LQR design
Q = 1*eye(nx);
R = 1*eye(nu);
[K,P,~] = dlqr(A,B,Q,R); % In MATLAB documentation, u(k)=-Kx(k)
K_lqr = -K; % In this lecture, we use the form u(k)=+Kx(k)

% (d) Closed-loop trajectory
k_sim = 100;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
for k = 1:k_sim
    uk(:,k) = K_lqr*xk(:,k);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

figure
subplot(1,3,1)
plot(0:k_sim,xk)
title('No saturation');
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');

% (e) LQR Stability
eig(A+B*K_lqr) % check if inside the unit circle

% (f) Input saturation
k_sim = 100;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
for k = 1:k_sim
    uk(:,k) = K_lqr*xk(:,k);
    uk(:,k) = min(max(uk(:,k),-1*ones(nu,1)),1*ones(nu,1)); % saturation -[1,1]
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

subplot(1,3,2)
plot(0:k_sim,xk)
title('Saturation [-1,1]');
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');

k_sim = 500;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
for k = 1:k_sim
    uk(:,k) = K_lqr*xk(:,k);
    uk(:,k) = min(max(uk(:,k),-0.1*ones(nu,1)),0.1*ones(nu,1)); % saturation -[0.1,0.1]
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

subplot(1,3,3)
plot(0:k_sim,xk)
title('Saturation [-0.1,0.1]');
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');

%% 3.2 Prediction model and cost function
N = 100;
% (a) Prediction model
% See ABN2PhiGamma.m
[Phi, Gamma] = ABN2PhiGamma(A,B,N);

% (b) Cost function
% See QRPN2PsiOmega.m
[Psi, Omega] = QRPN2PsiOmega(Q,R,P,N);
G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega*Phi;

%% 3.3 Unconstrained MPC
% (a) Unconstrained MPC feedback gain
K_mpc = -[eye(nu) zeros(nu,N-1)]*G^-1*F;

% (b) Closed-loop trajectory
k_sim = 100;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
for k = 1:k_sim
    uk(:,k) = K_mpc*xk(:,k);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

figure
plot(0:k_sim,xk)
title('Unconstrained MPC');
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');

% (c) Unconstrained MPC Stability
eig(A+B*K_mpc) % check if inside the unit circle

% (d) Compare LQR and unconstrained MPC