% 5LMB0 Model predictive control
% Instruction 4
clear all
close all

%% Initilization: Select one model
%First-order model

    Ac = 0;
    Bc = 1;
    Cc = 1;
    Dc = 0;
    x0 = 1; % Initial point
    Q = 10;
    R = 100;

% % Second-order model
%     Ac = [0 1;0 0];
%     Bc = [0;-1];
%     Cc = eye(2);
%     Dc = [0;0];
%     x0 = [3;1]; % Initial point
%     Q = [10,0;0,100];
%     R = 10;
%     x0 = [3;1]; % Initial point

% Discretization
Ts = 0.1;
ss_c = ss(Ac,Bc,Cc,Dc);
ss_d = c2d(ss_c,Ts,'zoh');
A = ss_d.A;
B = ss_d.B;
C = ss_d.C;
D = ss_d.D;
[nx,nu] = size(B);

%% 2.1 Terminal cost
% (a) LQR
[K,Pa,~] = dlqr(A,B,Q,R); % In MATLAB documentation, u(k)=-Kx(k)
Ka = -K; % In this lecture, we use the form u(k)=+Kx(k)

A_cl_LQR = A+B*Ka;

% (b) discrete Lyapunov equation
Kb = Ka;
Z = eye(nx);
% Z = Q+Kb'*R*Kb;
Pb = dlyap((A+B*Kb)',Z);

try chol(Pb)
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end

% (c) LMI
yalmip('clear')
O = sdpvar(nx,nx);
Y = sdpvar(nu,nx);
Con1 = [O, (A*O+B*Y)', O, Y';
    (A*O+B*Y), O, zeros(nx,nx), zeros(nx,nu);
    O, zeros(nx,nx), Q^-1, zeros(nx,nu);
    Y, zeros(nu,nx), zeros(nu,nx), R^-1]>=0;
Con2 = O>=1e-9;
constraints = Con1 + Con2;
diagnostics = optimize(constraints);
if diagnostics.problem == 0
    disp('Solver thinks it is feasible')
elseif diagnostics.problem == 1
    disp('Solver thinks it is infeasible')
else
    disp('Something else happened')
end
Pc = value(O)^-1;
Kc = value(Y)*value(O)^-1;

%% (d) Closed-loop trajectory
% In addition, you can check if J is strictly decreasing
close all
P = Pc; % Pa, Pb, Pc
N = 100;
[Phi, Gamma] = ABN2PhiGamma(A,B,N);
[Psi, Omega] = QRPN2PsiOmega(Q,R,P,N);
G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega*Phi;
K_mpc = -[eye(nu) zeros(nu,(N-1)*nu)]*G^-1*F;

k_sim = 500;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
for k = 1:k_sim
    uk(:,k) = K_mpc*xk(:,k);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

figure
stairs(0:k_sim,xk')
title('Unconstrained MPC');
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');
