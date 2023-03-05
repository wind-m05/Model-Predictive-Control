clear all, close all, clc

%% Simulation parameters
t = 300; 
Ts = 0.1;
N = 4; 
x0 = [10;-2.2541;-2.23*10^(-17);-0.2912];
u0 = [0.3792;
      0.0203];
% Reference
ref = zeros(2,length(t));
ref(1,1:200) = 10;
ref(1,201:end) = 8;
ref(2,1:100) = 0;
ref(2,101:end) = 5;
% Plotting
% stairs(t,ref(1,:))
% hold on
% stairs(t,ref(2,:))
% grid on
% xlabel('Time step [k]')
% ylabel('References')
% legend('speed v','climb rate h')

%% Predicted states bounds
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
Q = 10*eye(2);
R = 0.1*eye(nu);
T = diag(ones(N*nu,1));
n = size(T,1);
T(2:n+1:end) = -1;
% [x0,constr,Q,R] = defaulttest(nx,nu);

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
G = 2*(T'*psi*T+gamma'*C_bar'*omega*C_bar*gamma); % From Jelle
F = 2*gamma'*C_bar'*omega*C_bar*phi;
[W,L,c,S] = getWLcS(constr,N,B,gamma,phi,T);


v = [u0;zeros(2*(N-1),1)];
c = c+S*v; % update c to new constraints
% Update cost function here too !

%quadprog
xk = [x0 zeros(nx,t)];
uk = zeros(nu,t);
tk = zeros(1,t);
for k = 1:t
    tic;
    [Uk,fval,exitflag] = quadprog(G,F*xk(:,k),L,c+W*xk(:,k),[],[],[],[],[],[]);
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
 

% % 2. mpcActiveSetSolver
% xk = [x0 zeros(nx,t)];
% uk = zeros(nu,t);
% tk = zeros(1,t);
% Aeq = zeros(0,length(F));
% beq = zeros(0,1);
% [L_G,p] = chol((G+G')/2,'lower');
% Linv = linsolve(L_G,eye(size(L_G)),struct('LT',true));
% opt = mpcActiveSetOptions;
% opt.UseHessianAsInput = false;
% iA0 = false(size(c));
% for k = 1:t
%     tic
%     [Uk,exitflag,iA0,lambda] = mpcActiveSetSolver(Linv,F*xk(:,k),L,c+W*xk(:,k),Aeq,beq,iA0,opt);
%     tk(k) = toc;
%     uk(:,k) = Uk(1:nu);
%     xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
% end
% stairs(0:t-1,tk)
% legend('quadprog','mpcActiveSetSolver')


