% 5LMB0 Model predictive control
% Instruction 3
clear all
close all

%% 2.1.1 Feasible set of state and input: First-order model
Ac = 0;
Bc = 1;
Cc = 1;
Dc = 0;
Ts = 0.1;
ss_c = ss(Ac,Bc,Cc,Dc);
ss_d = c2d(ss_c,Ts,'zoh');
A = ss_d.A;
B = ss_d.B;
C = ss_d.C;
D = ss_d.D;
[nx,nu] = size(B);

xmin = -1;
xmax = 1;
umin = -1;
umax = 1;

N = 2;
x0 = 0.9;

[Phi, Gamma] = ABN2PhiGamma(A,B,N);
[W, L, c] = getWLc(A,B,xmax,xmin,umax,umin,Gamma,Phi);

% (a) Feasible set of inputs - LU<=c+Wx0
Feasibleset_U = Polyhedron('A',L,'B',c+W*x0);
figure
plot(Feasibleset_U)
xlabel('$u_{0|0}$','Interpreter','latex');
ylabel('$u_{1|0}$','Interpreter','latex');
sgtitle('Feasible set of inputs')

% (b) Feasible set of states - [-W L][x0;Uk]<=c
Feasibleset_x0_U = Polyhedron('A',[-W L],'B',c);
Feasibleset_x0 = projection(Feasibleset_x0_U,1,'vrep');
figure
plot(Feasibleset_x0)
xlabel('$x_{0|0}$','Interpreter','latex')
sgtitle('Feasible set of states')

%% 2.1.2 Feasible set of state and input: Second-order model
Ac = [0 1;0 0];
Bc = [0;-1];
Cc = eye(2);
Dc = [0;0];
Ts = 0.1;
ss_c = ss(Ac,Bc,Cc,Dc);
ss_d = c2d(ss_c,Ts,'zoh');
A = ss_d.A;
B = ss_d.B;
C = ss_d.C;
D = ss_d.D;
[nx,nu] = size(B);

xmin = [-10;-10];
xmax = [10;10];
umin = -0.1;
umax = 0.1;

N = 2;
x0 = [9;5];

[Phi, Gamma] = ABN2PhiGamma(A,B,N);
[W, L, c] = getWLc(A,B,xmax,xmin,umax,umin,Gamma,Phi);

% (a) Feasible set of inputs - LU<=c+Wx0
Feasibleset_U = Polyhedron('A',L,'B',c+W*x0);
figure
plot(Feasibleset_U)
xlabel('$u_{0|0}$','Interpreter','latex');
ylabel('$u_{1|0}$','Interpreter','latex');
sgtitle('Feasible set of inputs')

% (b) Feasible set of states - [-W L][x0;Uk]<=c
Feasibleset_x0_U = Polyhedron('A',[-W L],'B',c);
Feasibleset_x0 = projection(Feasibleset_x0_U,1:2,'vrep');
figure
plot(Feasibleset_x0)
xlabel('$x_{0|0}$','Interpreter','latex')
sgtitle('Feasible set of states')

%% 2.2 Constrained MPC
% (a) Feasibility problem will be discussed in Lecture 5

% (b)(c) Computation time
clear all
close all
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
[W, L, c] = getWLc(A,B,xmax,xmin,umax,umin,Gamma,Phi);

% Closed-loop trajectory
% 1. quadprog
k_sim = 500;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
opt =  optimoptions('quadprog','Display','off');
warning('off','optim:quadprog:HessianNotSym')
for k = 1:k_sim
    tic;
    [Uk,~,exitflag] = quadprog(G,F*xk(:,k),L,c+W*xk(:,k),[],[],[],[],[],opt);
    tk(k) = toc;
    if exitflag ~= 1
        warning('exitflag quadprog = %d\n', exitflag)
        if exitflag == -2
            sprintf('Optimization problem is infeasible.')
        end
    end
    uk(:,k) = Uk(1:nu);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

figure
hold on
xlim([0 k_sim]);
xlabel('$k$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
sgtitle('Computation time')
stairs(0:k_sim-1,tk)

% 2. mpcActiveSetSolver
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
Aeq = zeros(0,length(F));
beq = zeros(0,1);
[L_G,p] = chol((G+G')/2,'lower');
Linv = linsolve(L_G,eye(size(L_G)),struct('LT',true));
opt = mpcActiveSetOptions;
opt.UseHessianAsInput = false;
iA0 = false(size(c));
for k = 1:k_sim
    tic
    [Uk,exitflag,iA0,lambda] = mpcActiveSetSolver(Linv,F*xk(:,k),L,c+W*xk(:,k),Aeq,beq,iA0,opt);
    tk(k) = toc;
    uk(:,k) = Uk(1:nu);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end
stairs(0:k_sim-1,tk)
legend('quadprog','mpcActiveSetSolver')

%% 2.3 Explicit MPC
clear all
close all
% Second-order model
    Ac = [0 1;0 0];
    Bc = [0;-1];
    Cc = eye(2);
    Dc = [0;0];
    xmin = [-10;-10];
    xmax = [10;10];
    umin = -1;
    umax = 1;
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
N = 25;
[Phi, Gamma] = ABN2PhiGamma(A,B,N);
[Psi, Omega] = QRPN2PsiOmega(Q,R,P,N);
G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega*Phi;
[W, L, c] = getWLc(A,B,xmax,xmin,umax,umin,Gamma,Phi);

% (a) MPC Formulation using MPT3
model = LTISystem('A',A,'B',B);
model.x.min = xmin;
model.x.max = xmax;
% or use the following to set state constriant
% model.x.with('setConstraint');
% model.x.setConstraint = Polyhedron('lb', xmin, 'ub', xmax);
model.u.min = umin;
model.u.max = umax;
model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);
mpc = MPCController(model,N);
% to explicit
expmpc = mpc.toExplicit;
% expmpc.partition.plot()

% (b)(c)(d) Comparison with implicit MPC:
% 1. Explicit MPC
k_sim = 300;
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
for k = 1:k_sim
    tic;
    uk(:,k) = expmpc.evaluate(xk(:,k));
    tk(k) = toc;
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
sgtitle('Explicit MPC')

f_t = figure;
hold on
xlim([0 k_sim]);
xlabel('$k$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
sgtitle('Computation time')
stairs(0:k_sim-1,tk)

% 2. quadprog
% Compact formulation
xk = [x0 zeros(nx,k_sim)];
uk = zeros(nu,k_sim);
tk = zeros(1,k_sim);
opt =  optimoptions('quadprog','Display','off');
warning('off','optim:quadprog:HessianNotSym')
for k = 1:k_sim
    tic;
    [Uk,~,exitflag] = quadprog(G,F*xk(:,k),L,c+W*xk(:,k),[],[],[],[],[],opt);
    tk(k) = toc;
    if exitflag ~= 1
        warning('exitflag quadprog = %d\n', exitflag)
        if exitflag == -2
            sprintf('Optimization problem is infeasible.')
        end
    end
    uk(:,k) = Uk(1:nu);
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
sgtitle('Online constrained MPC')

figure(f_t);
hold on
xlim([0 k_sim]);
xlabel('$k$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
sgtitle('Computation time')
stairs(0:k_sim-1,tk)
legend('Explicit MPC','Online constrained MPC')

