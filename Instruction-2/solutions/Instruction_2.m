% 5LMB0 Model predictive control
% Instruction 2
clear all
close all

%% 2.1 Quadratic Programming
% (a)
% J = 1/2*2*u_0^2 + 5*u_0 + 25
%     s.t. [1;-1]*u_0<=[1;1] or -1<=u_0<=1
u0_1 = quadprog(2,5,[1;-1],[1;1]) %or u_0 = quadprog(2,5,[],[],-1,1)

% (b)
% J = 1/2*9/4*u_0^2 + 7/2*u_0 + 9
%     s.t. [1;-1]*u_0<=[1;1] or -1<=u_0<=1
u0_2 = quadprog(9/4,7/2,[1;-1],[1;1]) %or u_0 = quadprog(9/4,7/2,[],[],-1,1)

clear u0_1 u0_2
clc

%% 2.2 Constrained MPC
% Initilization: Select one model
% First-order model
%     Ac = 0;
%     Bc = 1;
%     Cc = 1;
%     Dc = 0;
%     xmin = -1e9;
%     xmax = 1e9;
%     umin = -0.3;
%     umax = 0.3;
%     x0 = 1; % Initial point

% Second-order model
    Ac = [0 1;0 0];
    Bc = [0;-1];
    Cc = eye(2);
    Dc = [0;0];
    xmin = [-10;-10];
    xmax = [10;10];
    umin = -0.1;
    umax = 0.1;
    x0 = [3;1]; % Initial point

% Discretization
Ts = 0.1;
ss_c = ss(Ac,Bc,Cc,Dc);
ss_d = c2d(ss_c,Ts,'zoh');
A = ss_d.A;
B = ss_d.B;
C = ss_d.C;
D = ss_d.D;
[nx,nu] = size(B);

% Q,R,P design
Q = 1*eye(nx);
R = 1*eye(nu);
[K,P,~] = dlqr(A,B,Q,R); % In MATLAB documentation, u(k)=-Kx(k)
K_lqr = -K; % In this lecture, we use the form u(k)=+Kx(k)

% Compact formulation
N = 100;
[Phi, Gamma] = ABN2PhiGamma(A,B,N);
[Psi, Omega] = QRPN2PsiOmega(Q,R,P,N);
G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega*Phi;

% (a)(b) Constraints: Compact formulation
% See getWLc.m
[W, L, c] = getWLc(A,B,xmax,xmin,umax,umin,Gamma,Phi);

% (c) Closed-loop trajectory
k_sim = 500;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
Uk = zeros(nu*N,k_sim);
% opt =  optimoptions('quadprog','Display','off');
% warning('off','optim:quadprog:HessianNotSym')
for k = 1:k_sim
    [U,~,exitflag] = quadprog(G,F*xk(:,k),L,c+W*xk(:,k),[],[],[],[],[],[]);
    if exitflag ~= 1
        warning('exitflag quadprog = %d\n', exitflag)
        if exitflag == -2
            sprintf('Optimization problem is infeasible.')
        end
    end
    Uk(:,k) = U;
    uk(:,k) = U(1:nu);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

figure
subplot(1,2,1);
stairs(0:k_sim,xk')
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');
xlim([0 k_sim]);

subplot(1,2,2);
stairs(0:k_sim-1,uk)
xlabel('$k$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xlim([0 k_sim]);
ylim([umin-0.1 umax+0.1]);
sgtitle('Constrained MPC')

% (d) Comparison with open loop predicted trajecotry
% Follow the previous section
figure
hold on
for k = 1:k_sim
    stairs(k-1:k+N-2,Uk(:,k))
end
stairs(0:k_sim-1,uk,'k')
xlabel('$k$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xlim([0 k_sim]);
ylim([umin-0.1 umax+0.1]);
title('Open loop predicted trajecotry of input')

% (e) Comparison with saturated unconstrained MPC
K_mpc = -[eye(nu) zeros(nu,N-1)]*G^-1*F;
k_sim = 500;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
for k = 1:k_sim
    uk(:,k) = K_mpc*xk(:,k);
    uk(:,k) = min(max(uk(:,k),umin),umax); % saturation -[1,1]
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

figure
subplot(1,2,1);
stairs(0:k_sim,xk')
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');
xlim([0 k_sim]);

subplot(1,2,2);
stairs(0:k_sim-1,uk)
xlabel('$k$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xlim([0 k_sim]);
ylim([umin-0.1 umax+0.1]);
sgtitle('Saturated unconstrained MPC')