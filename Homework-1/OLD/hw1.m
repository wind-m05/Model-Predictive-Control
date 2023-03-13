% clear all, close all, clc
%% Simulation parameters
t = 300; 
Ts = 0.1;
N = 20; 
x0 = [10;-2.2541;-2.23*10^(-17);-0.2912];
u0 = [0.3792;
      0.0203];
y0 = [10;0];

ref = zeros(2,t+N);
ref(1,1:200) = 10;
ref(1,201:end) = 8;
ref(2,1:100) = 0;
ref(2,101:end) = 5;
ref_test = ref;
ref = reshape(ref,2*(t+N),1);
%%
% Predicted states bounds
constr.statelb = [-25;-25;-pi/20;-pi/2;-15]; % v,w,q,theta,h
constr.stateub =  [25;25;pi/20;pi/2;15]; % v,w,q,theta,h
constr.initialstatelb = -abs([x0; 0]); % v,w,q,theta,h
constr.initialstateub = abs([x0; 0]);
constr.terminalstatelb = [-25;-25;-pi/20;-pi/2;-15]; % Change if needed after monday
constr.terminalstateub = [25;25;pi/20;pi/2;15];
constr.inputlb = [-20;-20];
constr.inputub =  [20;20];
constr.deltainputlb = [-5;-5];
constr.deltainputub =  [5;5];

% Model
[A,B,C,D,sys] = modelselect('aircraft','discrete',Ts);
[nx,nu] = size(B);
ny = size(C,1);
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

[phi,gamma] = phi_gam_predict(A,B,N);
[omega,psi] = get_omega_psi(Q,Q,R,N); % If we have a P change second Q into P
Cc = repmat({C}, 1, N);  
C_bar = blkdiag(Cc{:});
[W,L,c,S] = getWLcS(constr,N,B,gamma,phi,T);


%% Unconstrained
% constraint = 'unconstrained'; 
% uk = [u0,zeros(nu,t)];
% xk = [x0 A*x0+B*u0 zeros(nx,t)];
% yk = [y0 C*xk(:,2) zeros(ny,t)];
% 
% v = [u0 ; zeros(2*(N-1),1)];
% G = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
% M = zeros(size(B,2),2*N);
% M(1:size(B,2),1:size(B,2)) = eye(size(B,2));
% i = 0;
% for k = 2:(t+1)
%     i = i+2; 
%     Rk = ref((i+1):(i+ny*N));
%     v = [uk(:,k-1); zeros(2*(N-1),1)];
%     uk(:,k) = -M*inv(G)*2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,k)-gamma'*C_bar'*omega*Rk-T'*psi*v); 
%     xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
%     yk(:,k) = C*xk(:,k);
% end

%% Constrained quadprog
% constraint = 'constrained quadprog'; 
% uk = [u0,zeros(nu,t)];
% xk = [x0 A*x0+B*u0 zeros(nx,t-1)];
% yk = [y0 C*xk(:,2) zeros(ny,t-1)];
% H = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
% i = 0;
% for k = 2:t
%     i = i+2; 
%     Rk = ref((i+1):(i+ny*N));
%     v = [uk(:,k-1) ; zeros(2*(N-1),1)];
%     f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,k)-gamma'*C_bar'*omega*Rk-T'*psi*v);
%     d = c + S*v;
%     [Uk,fval,exitflag] = quadprog(H,f,L,d+W*xk(:,k),[],[],[],[],[],[]);
%     if exitflag ~= 1
%         warning('exitflag quadprog = %d\n', exitflag)
%         if exitflag == -2
%             sprintf('Optimization problem is infeasible.')
%         end
%     end
%     uk(:,k) = Uk(1:nu);
%     xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
%     yk(:,k+1) = C*xk(:,k+1);
% end



% mpcActiveSetSolver
constraint = 'constrained mpcsetsolver'; 
uk = [u0,zeros(nu,t)];
xk = [x0 A*x0+B*u0 zeros(nx,t-1)];
yk = [y0 C*xk(:,2) zeros(ny,t-1)];
v = [u0 ; zeros(2*(N-1),1)];
H = 2*(gamma'*C_bar'*omega*C_bar*gamma+T'*psi*T);
f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,1)-gamma'*C_bar'*omega*ref(1:2*N)-2*T'*psi*v);
G = H;
Aeq = zeros(0,length(f));
beq = zeros(0,1);
[L_G,p] = chol((G+G')/2,'lower');
Linv = linsolve(L_G,eye(size(L_G)),struct('LT',true));
opt = mpcActiveSetOptions;
opt.UseHessianAsInput = false;
iA0 = false(size(c));
i = 0;
for k = 2:t
    i = i+2; 
    Rk = ref((i+1):(i+ny*N));
    f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,k)-gamma'*C_bar'*omega*Rk-T'*psi*v);
    d = c + S*v;
    v = [uk(:,k-1) ; zeros(2*(N-1),1)];
    [Uk,exitflag,iA0,lambda] = mpcActiveSetSolver(Linv,f,L,d+W*xk(:,k),Aeq,beq,iA0,opt);
    uk(:,k) = Uk(1:nu);
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
    yk(:,k+1) = C*xk(:,k+1);
    d = c + S*v;
    f = 2*(gamma'*C_bar'*omega*C_bar*phi*xk(:,1)-gamma'*C_bar'*omega*ref(1:2*N)-2*T'*psi*v);
end

ref = reshape(ref,[2,t+N]);
ref1 = ref(1,:);
ref2 = ref(2,:);
font = 18;
thickness = 2;
plottingFunction_michiel(constraint,font,thickness,t,xk,uk,yk,ref1,ref2);
t_mic = t;
xk_mic = xk;
uk_mic = uk;
yk_mic = yk;
ref1_mic = ref1;
ref2_mic = ref2;
save('michiel.mat','t_mic','xk_mic','uk_mic','yk_mic','ref1_mic','ref2_mic')

