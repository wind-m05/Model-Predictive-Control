% Model predictive control Homework assignment 1
% Jelle Cruijsen & Michiel Wind
close all;
clear all;
clc

t = 300; 
Ts = 0.1;
N = 20;

%% Question 4.1 (discrete-time state-space prediction model)
Ac = [-0.003 0.039 0 -0.322; -0.065 -0.319 7.74 0; 0.02 -0.101 -0.429 0; 0 0 1 0];
Bc = [0.01 1; -0.18 -0.04; -1.16 0.598; 0 0];
Cc = [1 0 0 0; 0 -1 0 7.74];
Dc = 0;

ss_c = ss(Ac,Bc,Cc,Dc);
ss_d = c2d(ss_c,Ts,'zoh');
A = ss_d.A;
B = ss_d.B;
C = ss_d.C;
D = ss_d.D;
[n,m] = size(B);
p = height(C);

% Initial conditions
h0 = 0;
v0 = 10;
x0 = [v0; -2.2541; -2.23e-17; -0.2912];
u0 = [0.3792; 0.0203];
y0 = [v0; h0];

ref = zeros(2*(t+N),1);
k = 0;
for i = 1:2:length(ref)
    if(k<100)
        ref(i) = 10;
        ref(i+1) = 0;
    elseif(k<200)
        ref(i) = 10;
        ref(i+1) = 5;
    else
        ref(i) = 8;
        ref(i+1) = 5;
    end
    k = k+1;
end
ref1 = zeros((t+N),1);
k = 0;
for i = 1:2:length(ref)
    k = k+1;
    ref1(k) = ref(i);
end
ref2 = zeros((t+N),1);
k = 0;
for i = 1:2:length(ref)
    k = k+1;
    ref2(k) = ref(i+1);
end

%% Constraints and matrices
umin = [-20; -20];
umax = [20; 20];
xmin = [-25; -25; -pi/20; -pi/2];
xmax = [25; 25; pi/20; pi/2];
ymin = -15;
ymax = 15;
dumin = [-15; -15];
dumax = [15; 15];

Q = 10*eye(p);
R = 0.1*eye(m);

[phi, gamma] = predictionModel(A,B,N,n,m);
omega = kron(eye(N),Q);       % Kronecker product
psi = kron(eye(N),R);
C_bar = kron(eye(N),C);

T = diag(ones(N*m,1));
nT = size(T,1);
T(2:nT+1:end) = -1;

%% Question 4.2 (design a model predictive controller)
k_sim = t;
% unconstrained
uk = [u0 zeros(m,t)];
xk = [x0 A*x0+B*u0 zeros(n,t)];
yk = [y0 C*xk(:,2) zeros(p,t-1)];
G = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
i = 0;
for k = 2:(k_sim+1)               %time starts at k = 0, xk is known for k = 1
    i = i+2;                      %used for indexing due to reference being twice as long;
    Rk = ref((i+1):(i+p*N));
    v = [uk(:,k-1); zeros(2*(N-1),1)];
    Uk = -2*G^-1*(gamma'*C_bar'*omega*C_bar*phi*xk(:,k)-gamma'*C_bar'*omega*Rk-T'*psi*v);
    uk(:,k) = Uk(1:m);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
    yk(:,k) = C*xk(:,k);
end
x_unconstrained = xk;
%% constrained
[Ccal, Dcal, Ecal, Mcal,Ebar] = caligraphicMatricesExtended(umin,umax,xmin,xmax,ymin,ymax,dumin,dumax,N,p,n,m);  
uk = [u0 zeros(m,t)];
xk = [x0 A*x0+B*u0 zeros(n,t)];
yk = [y0 C*xk(:,2) zeros(p,t-1)];
H = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
L = Mcal*gamma + Ecal;  %no delta u and y constraints
L = Mcal*gamma + Ecal+Ebar*T;
W = -Dcal-Mcal*phi;
i = 0;
for k = 2:(k_sim+1)               %time starts at k = 0, xk is known for k = 1
    i = i+2;                      %used for indexing due to reference being twice as long;
    Rk = ref((i+1):(i+p*N));
    v = [uk(:,k-1); zeros(2*(N-1),1)];
    f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,k)-gamma'*C_bar'*omega*Rk-T'*psi*v);
%     [Uk,fval,exitflag] = quadprog(H,f,L,Ccal+W*xk(:,k),[],[],[],[],[],[]);  %no delta u and y constraints
    [Uk,fval,exitflag] = quadprog(H,f,L,Ccal+W*xk(:,k)+Ebar*v,[],[],[],[],[],[]);
    if exitflag ~= 1
        warning('exitflag quadprog = %d\n', exitflag)
        if exitflag == -2
            sprintf('Optimization problem is infeasible.')
            break;                %optimization failed, break loop then plot results
        end
    end
    uk(:,k) = Uk(1:m);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
    yk(:,k) = C*xk(:,k);
end

% %% 2. mpcActiveSetSolver
% [Ccal, Dcal, Ecal, Mcal] = caligraphicMatrices(umin,umax,xmin,xmax,N,n,m);  %no delta u and y constraints
% [Ccal, Dcal, Ecal, Mcal,Ebar] = caligraphicMatricesExtended(umin,umax,xmin,xmax,ymin,ymax,dumin,dumax,N,p,n,m);  
% uk = [u0 zeros(m,t)];
% xk_as = [x0 A*x0+B*u0 zeros(n,t)];
% yk = [y0 C*xk_as(:,2) zeros(p,t-1)];
% v = [uk(:,1); zeros(2*(N-1),1)];
% Rk = ref((2+1):(2+p*N));
% G = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
% f = 2*(gamma'*C_bar'*omega*C_bar*phi*x0-gamma'*C_bar'*omega*Rk-T'*psi*v);
% Aeq = zeros(0,length(f));
% beq = zeros(0,1);
% [L_G,p_chol] = chol((G+G')/2,'lower');
% Linv = linsolve(L_G,eye(size(L_G)),struct('LT',true));
% opt = mpcActiveSetOptions;
% opt.UseHessianAsInput = false;
% iA0 = false(size(Ccal));
% i=0;
% L = Mcal*gamma + Ecal+Ebar*T;
% W = -Dcal-Mcal*phi;
% for k = 2:(k_sim+1)
%     i = i+2;                      %used for indexing due to reference being twice as long;
%     Rk = ref((i+1):(i+p*N));
%     v = [uk(:,k-1); zeros(2*(N-1),1)];
%     f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk_as(:,k)-gamma'*C_bar'*omega*Rk-T'*psi*v);
%     [Uk,exitflag,iA0,lambda] = mpcActiveSetSolver(Linv,f,L,Ccal+W*xk_as(:,k)+Ebar*v,Aeq,beq,iA0,opt);
% %     [Uk,fval,exitflag] = quadprog(H,f,L,Ccal+W*xk(:,k)+Ebar*v,[],[],[],[],[],[]);
%     uk(:,k) = Uk(1:m);
%     xk_as(:,k+1) = A*xk(:,k)+B*uk(:,k);
%     yk(:,k) = C*xk(:,k);
% end


%% plotting
figure()
subplot(1,2,1)
stairs(0:t,yk(1,:))
xlabel('$k$','Interpreter','latex');
ylabel('$v [ft/sec]$','Interpreter','latex');

subplot(1,2,2)
stairs(0:t,yk(2,:))
xlabel('$k$','Interpreter','latex');
ylabel('$h [ft/sec]$','Interpreter','latex');
sgtitle('Outputs')

figure()
subplot(2,2,1)
stairs(0:t,xk(1,1:length(xk)-1))
xlabel('$k$','Interpreter','latex');
ylabel('$v [ft/sec]$','Interpreter','latex');

subplot(2,2,2)
stairs(0:t,xk(2,1:length(xk)-1))
xlabel('$k$','Interpreter','latex');
ylabel('$w [ft/sec]$','Interpreter','latex');

subplot(2,2,3)
stairs(0:t,xk(3,1:length(xk)-1))
xlabel('$k$','Interpreter','latex');
ylabel('$q [ft/sec]$','Interpreter','latex');

subplot(2,2,4)
stairs(0:t,xk(4,1:length(xk)-1))
xlabel('$k$','Interpreter','latex');
ylabel('$\theta [rad/sec]$','Interpreter','latex');
sgtitle('States')


figure()
subplot(1,2,1)
stairs(0:t,uk(1,:))
xlabel('$k$','Interpreter','latex');
ylabel('$e$','Interpreter','latex');

subplot(1,2,2)
stairs(0:t,uk(2,:))
xlabel('$k$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex');
sgtitle('Inputs')

figure()
subplot(1,2,1)
stairs(0:t-1,ref1(1:t))
xlabel('$k$','Interpreter','latex');
ylabel('$v ref$','Interpreter','latex');
ylim([min(ref1)-1 max(ref1)+1]);

subplot(1,2,2)
stairs(0:t-1,ref2(1:t))
xlabel('$k$','Interpreter','latex');
ylabel('$h ref$','Interpreter','latex');
sgtitle('reference');
ylim([min(ref2)-1 max(ref2)+1]);

%% Comparison quadprog/active set solver

% figure()
% subplot(2,2,1)
% stairs(0:t,xk(1,1:length(xk)-1))
% hold on
% stairs(0:t,xk_as(1,1:length(xk_as)-1))
% xlabel('$k$','Interpreter','latex');
% ylabel('$v [ft/sec]$','Interpreter','latex');
% 
% subplot(2,2,2)
% stairs(0:t,xk(2,1:length(xk)-1))
% hold on
% stairs(0:t,xk_as(2,1:length(xk_as)-1))
% xlabel('$k$','Interpreter','latex');
% ylabel('$w [ft/sec]$','Interpreter','latex');
% 
% subplot(2,2,3)
% stairs(0:t,xk(3,1:length(xk)-1))
% hold on
% stairs(0:t,xk_as(3,1:length(xk_as)-1))
% xlabel('$k$','Interpreter','latex');
% ylabel('$q [ft/sec]$','Interpreter','latex');
% 
% subplot(2,2,4)
% stairs(0:t,xk(4,1:length(xk)-1))
% hold on
% stairs(0:t,xk_as(4,1:length(xk_as)-1))
% xlabel('$k$','Interpreter','latex');
% ylabel('$\theta [rad/sec]$','Interpreter','latex');
% sgtitle('States')

%% Comparison constrained/unconstrained

figure()
subplot(2,2,1)
stairs(0:t,x_unconstrained(1,1:length(x_unconstrained)-1))
hold on
stairs(0:t,xk(1,1:length(xk)-1))
xlabel('$k$','Interpreter','latex');
ylabel('$v [ft/sec]$','Interpreter','latex');
legend('unconstrained','constrained')
subplot(2,2,2)
stairs(0:t,x_unconstrained(2,1:length(x_unconstrained)-1))
hold on
stairs(0:t,xk(2,1:length(xk)-1))
xlabel('$k$','Interpreter','latex');
ylabel('$w [ft/sec]$','Interpreter','latex');

subplot(2,2,3)
stairs(0:t,x_unconstrained(3,1:length(x_unconstrained)-1))
hold on
stairs(0:t,xk(3,1:length(xk)-1))
xlabel('$k$','Interpreter','latex');
ylabel('$q [ft/sec]$','Interpreter','latex');

subplot(2,2,4)
stairs(0:t,x_unconstrained(4,1:length(x_unconstrained)-1))
hold on
stairs(0:t,xk(4,1:length(xk)-1))
xlabel('$k$','Interpreter','latex');
ylabel('$\theta [rad/sec]$','Interpreter','latex');
sgtitle('States')
