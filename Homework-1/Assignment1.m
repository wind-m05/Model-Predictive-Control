% Model predictive control Homework assingment 1
% Jelle Cruijsen & Michiel Wind

close all; clear all; clc

% Simulation variables
t = 300; 
Ts = 0.1;
N = 20;

% Model definition
[A,B,C,D,sys] = modelselect('aircraft','discrete',Ts);
[n,m] = size(B);
p = height(C);

% Initial conditions
h0 = 0;
v0 = 10;
x0 = [v0; -2.2541; -2.23e-17; -0.2912];
u0 = [0.3792; 0.0203];
y0 = [v0; h0];

% Reference definition
ref = zeros(2,t+N);
ref(1,1:200) = 10;
ref(1,201:end) = 8;
ref(2,1:100) = 0;
ref(2,101:end) = 5;
ref = reshape(ref,2*(t+N),1);

% Constraints
umin = [-20; -20];
umax = [20; 20];
xmin = [-25; -25; -pi/20; -pi/2];
xmax = [25; 25; pi/20; pi/2];
ymin = -15;
ymax = 15;
dumin = [-15; -15];
dumax = [15; 15];

% Performance
Q = 10*eye(p);
R = 0.1*eye(m);

[phi, gamma] = predictionModel(A,B,N,n,m);
omega = kron(eye(N),Q);       % Kronecker product
psi = kron(eye(N),R);
C_bar = kron(eye(N),C);

T = diag(ones(N*m,1)); % Delta matrix to handle delta constraints
nT = size(T,1);
T(2:nT+1:end) = -1;

%% unconstrained MPC
% constraint = 'unconstrained';           %changes plot title
% uk = [u0 zeros(m,t)];
% xk = [x0 A*x0+B*u0 zeros(n,t-1)];
% yk = [y0 C*xk(:,2) zeros(p,t-1)];
% G = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
% i = 0;
% for k = 2:t              % time starts at k = 0, xk is known for k = 1
%     i = i+2;             % used for indexing due to reference being twice as long;
%     Rk = ref((i+1):(i+p*N));
%     v = [uk(:,k-1); zeros(2*(N-1),1)];
%     Uk = -2*G^-1*(gamma'*C_bar'*omega*C_bar*phi*xk(:,k)-gamma'*C_bar'*omega*Rk-T'*psi*v);
%     uk(:,k) = Uk(1:m);
%     xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
%     yk(:,k+1) = C*xk(:,k+1);
% end

%% constrained MPC Quadprog
constraint = 'constrained';    
[Ccal, Dcal, Ecal, Mcal,Ebar] = caligraphicMatricesExtended(umin,umax,xmin,xmax,ymin,ymax,dumin,dumax,N,p,n,m);
% [D,M,E,S,c] = DMESc_constraints(constr,N,B);
% Mcal = M;
uk = [u0 zeros(m,t)];
xk = [x0 A*x0+B*u0 zeros(n,t-1)];
yk = [y0 C*xk(:,2) zeros(p,t-1)];
H = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
L = Mcal*gamma + Ecal+Ebar*T;
W = -Dcal-Mcal*phi;
i = 0;
for k = 2:t              %time starts at k = 0, xk is known for k = 1
    i = i+2;                      %used for indexing due to reference being twice as long;
    Rk = ref((i+1):(i+p*N));
    v = [uk(:,k-1); zeros(2*(N-1),1)];
    f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,k)-gamma'*C_bar'*omega*Rk-T'*psi*v);
    [Uk,fval,exitflag] = quadprog(H,f,L,Ccal+W*xk(:,k)+Ebar*v,[],[],[],[],[],[]);
    if exitflag ~= 1
        warning('exitflag quadprog = %d\n', exitflag)
        if exitflag == -2
            sprintf('Optimization problem is infeasible.')
            break;  %optimization failed, break loop then plot results
        end
    end
    uk(:,k) = Uk(1:m);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
    yk(:,k+1) = C*xk(:,k+1);
end

%% Constrained MPC active set solver
% nu = m;
% nx = n;
% ny = p;
% [Ccal, Dcal, Ecal, Mcal,Ebar] = caligraphicMatricesExtended(umin,umax,xmin,xmax,ymin,ymax,dumin,dumax,N,p,n,m); 
% c = Ccal;
% S = Ecal;
% L = Mcal*gamma + Ecal+Ebar*T;
% W = -Dcal-Mcal*phi;
% constraint = 'constrained mpcsetsolver'; 
% uk = [u0,zeros(nu,t)];
% xk = [x0 A*x0+B*u0 zeros(nx,t-1)];
% yk = [y0 C*xk(:,2) zeros(ny,t-1)];
% v = [u0 ; zeros(2*(N-1),1)];
% H = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
% f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,1)-gamma'*C_bar'*omega*ref(1:2*N)-2*T'*psi*v);
% G = H;
% Aeq = zeros(0,length(f));
% beq = zeros(0,1);
% [L_G,p] = chol((G+G')/2,'lower');
% Linv = linsolve(L_G,eye(size(L_G)),struct('LT',true));
% opt = mpcActiveSetOptions;
% opt.UseHessianAsInput = false;
% iA0 = false(size(c));
% i = 0;
% for k = 2:t
%     i = i+2; 
%     Rk = ref((i+1):(i+ny*N));
%     f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,k)-gamma'*C_bar'*omega*Rk-T'*psi*v);
%     d = c + S*v;
%     v = [uk(:,k-1) ; zeros(2*(N-1),1)];
%     [Uk,exitflag,iA0,lambda] = mpcActiveSetSolver(Linv,f,L,d+W*xk(:,k),Aeq,beq,iA0,opt);
%     uk(:,k) = Uk(1:nu);
%     xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
%     yk(:,k+1) = C*xk(:,k+1);
%     d = c + S*v;
%     f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,1)-gamma'*C_bar'*omega*ref(1:2*N)-2*T'*psi*v);
% end

%% plotting
ref = reshape(ref,[2,t+N]);
ref1 = ref(1,:);
ref2 = ref(2,:);
font = 18;
thickness = 2;
plottingFunction(constraint,font,thickness,t,xk,uk,yk,ref1,ref2);
