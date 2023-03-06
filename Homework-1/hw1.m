clear all, close all, clc
%% Simulation parameters
t = 300; 
Ts = 0.1;
N = 20; 
x0 = [10;-2.2541;-2.23*10^(-17);-0.2912];
u0 = [0.3792;
      0.0203];
y0 = [10;0];

% Reference

ref = zeros(2,t+N);
ref(1,1:200) = 10;
ref(1,201:end) = 8;
ref(2,1:100) = 0;
ref(2,101:end) = 5;
ref = reshape(ref,2*(t+N),1);

% Predicted states bounds
constr.statelb = [-25;-25;-pi/20;-pi/2;-15]; % v,w,q,theta,h
constr.stateub =  [25;25;pi/20;pi/2;15]; % v,w,q,theta,h
constr.initialstatelb = -abs([x0; 0]); % v,w,q,theta,h
constr.initialstateub = abs([x0; 0]);
constr.terminalstatelb = [-25;-25;-pi/20;-pi/2;-15]; % Change if needed after monday
constr.terminalstateub = [25;25;pi/20;pi/2;15];
constr.inputlb = [-20;-20];
constr.inputub =  [20;20];
constr.deltainputlb = [-15;-15];
constr.deltainputub =  [15;15];

% Model
[A,B,C,D,sys] = modelselect('aircraft','discrete',Ts);
[nx,nu] = size(B);
ny = 2;
Q = 10*eye(ny);
R = 0.1*eye(nu);
T = diag(ones(N*nu,1));
n = size(T,1);
T(2:n+1:end) = -1;

%% Checks
ctrb_M = ctrb(A,B);
if length(x0) ~= length(A)
error('length of initial conditions must be the same as A')
elseif rank(ctrb_M) ~= length(A)
error('Problem is uncontrollable') 
elseif length(R) ~= nu
error('R matrix must be same size as number of inputs')
end

% [Klqr,P,~] = dlqr(A,B,Q,R); % CHECK with lecturer if we need P in this
% exercise...
[phi,gamma] = phi_gam_predict(A,B,N);
[omega,psi] = get_omega_psi(Q,Q,R,N); % If we have a P change second Q into P
Cc = repmat({C}, 1, N);  
C_bar = blkdiag(Cc{:});
% G = 2*(T'*psi*T+gamma'*C_bar'*omega*C_bar*gamma); % From Jelle
% F = 2*gamma'*C_bar'*omega*C_bar*phi;
[W,L,c,S] = getWLcS(constr,N,B,gamma,phi,T);


%% Unconstrained
xk = [x0 zeros(nx,t)];
yk = [y0 zeros(ny,t)];
uk = [zeros(size(u0)),u0,zeros(nu,t-1)];
xk(:,2) = A*xk(:,1)+B*uk(:,1); % Calculate x1 (because we know u0)
yk(:,2) = C*xk(:,1);
v = [u0 ; zeros(2*(N-1),1)];
G = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
F = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,1)-gamma'*C_bar'*omega*ref(1:2*N)-T'*psi*v);
M = zeros(size(B,2),2*N);
M(1:size(B,2),1:size(B,2)) = eye(size(B,2 ));
for k = 2:t
    uk(:,k) = -M*inv(G)*(gamma'*C_bar'*omega*C_bar*phi*xk(:,k)-gamma'*C_bar'*omega*ref(k:(k-1)+(2*N))-T'*psi*v); 
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
    yk(:,k+1) = C*xk(:,k);
    v = [uk(:,k) ; zeros(2*(N-1),1)];
end

% %% quadprog
% xk = [x0 zeros(nx,t)];
% yk = [y0 zeros(ny,t)];
% uk = [zeros(size(u0)),u0,zeros(nu,t-1)];
% xk(:,2) = A*xk(:,1)+B*uk(:,1); % Calculate x1 (because we know u0)
% yk(:,2) = C*xk(:,1);
% v = [u0 ; zeros(2*(N-1),1)];
% H = (gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
% f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,1)-gamma'*C_bar'*omega*ref(1:2*N)-2*T'*psi*v);
% c = c+S*v; % update c to new constraints
% for k = 2:t
%     [Uk,fval,exitflag] = quadprog(H,f,L,c+W*xk(:,k),[],[],[],[],[],[]);
%     if exitflag ~= 1
%         warning('exitflag quadprog = %d\n', exitflag)
%         if exitflag == -2
%             sprintf('Optimization problem is infeasible.')
%         end
%     end
%     uk(:,k) = Uk(1:nu);
%     xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
%     yk(:,k+1) = C*xk(:,k);
%     v = [uk(:,k) ; zeros(2*(N-1),1)];
%     c = c+S*v; 
%     f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,1)-gamma'*C_bar'*omega*ref(k:(k-1)+(2*N))-2*T'*psi*v);
% end

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
subplot(1,2,1)
stairs(0:t,uk(1,:))
xlabel('$k$','Interpreter','latex');
ylabel('$e$','Interpreter','latex');

subplot(1,2,2)
stairs(0:t,uk(2,:))
xlabel('$k$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex');
sgtitle('Inputs')

ref = reshape(ref,[2,t+N]);
figure()
subplot(1,2,1)
stairs(0:t-1+N,ref(1,:)')
xlabel('$k$','Interpreter','latex');
ylabel('$v ref$','Interpreter','latex');

subplot(1,2,2)
stairs(0:t-1+N,ref(2,:)')
xlabel('$k$','Interpreter','latex');
ylabel('$h ref$','Interpreter','latex');
sgtitle('reference')


 
% %% mpcActiveSetSolver
% xk = [x0 zeros(nx,t)];
% yk = [y0 zeros(ny,t)];
% uk = [zeros(size(u0)),u0,zeros(nu,t-1)];
% xk(:,2) = A*xk(:,1)+B*uk(:,1); % Calculate x1 (because we know u0)
% yk(:,2) = C*xk(:,1);
% v = [u0 ; zeros(2*(N-1),1)];
% H = (gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
% f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,1)-gamma'*C_bar'*omega*ref(1:2*N)-2*T'*psi*v);
% c = c+S*v; % update c to new constraints
% G = H;
% c = c+S*v; % update c to new constraints
% Aeq = zeros(0,length(f));
% beq = zeros(0,1);
% [L_G,p] = chol((G+G')/2,'lower');
% Linv = linsolve(L_G,eye(size(L_G)),struct('LT',true));
% opt = mpcActiveSetOptions;
% opt.UseHessianAsInput = false;
% iA0 = false(size(c));
% for k = 2:t
%     [Uk,exitflag,iA0,lambda] = mpcActiveSetSolver(Linv,f,L,c+W*xk(:,k),Aeq,beq,iA0,opt);
%     uk(:,k) = Uk(1:nu);
%     xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
%     yk(:,k+1) = C*xk(:,k);
%     v = [uk(:,k) ; zeros(2*(N-1),1)];
%     c = c+S*v; 
%     f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,1)-gamma'*C_bar'*omega*ref(1:2*N)-2*T'*psi*v);
% end
% 
% figure()
% subplot(1,2,1)
% stairs(0:t,yk(1,:))
% xlabel('$k$','Interpreter','latex');
% ylabel('$v [ft/sec]$','Interpreter','latex');
% 
% subplot(1,2,2)
% stairs(0:t,yk(2,:))
% xlabel('$k$','Interpreter','latex');
% ylabel('$h [ft/sec]$','Interpreter','latex');
% sgtitle('Outputs')
% 
% figure()
% subplot(1,2,1)
% stairs(0:t,uk(1,:))
% xlabel('$k$','Interpreter','latex');
% ylabel('$e$','Interpreter','latex');
% 
% subplot(1,2,2)
% stairs(0:t,uk(2,:))
% xlabel('$k$','Interpreter','latex');
% ylabel('$\tau$','Interpreter','latex');
% sgtitle('Inputs')


