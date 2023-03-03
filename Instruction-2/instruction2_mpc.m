clear all, close all, clc

% Simulation parameters
t = 500; 
Ts = 0.1;
N = 100; 
x0 = [3;1];

% Predicted states bounds
constr.statelb = [-10;-10]; 
constr.stateub =  [10;10];
constr.initialstatelb = -abs(x0);
constr.initialstateub = abs(x0);
constr.terminalstatelb = [-3;-3];
constr.terminalstateub = [3;3];
constr.inputlb = -0.1;
constr.inputub =  0.1;

% Model
[A,B,C,D,sys] = modelselect('car2d','discrete',Ts);
[nx,nu] = size(B);
Q = 1*eye(nx);
R = 1*eye(nu);
% [x0,constr,Q,R] = defaulttest(nx,nu);


% Checks
ctrb_M = ctrb(A,B);
if length(x0) ~= length(A)
error('length of initial conditions must be the same as A')
elseif rank(ctrb_M) ~= length(A)
error('Problem is uncontrollable') 
elseif length(Q) ~= nx
error('Q matrix must be same size as number of states') 
elseif length(R) ~= nu
error('R matrix must be same size as number of inputs')
elseif length(constr.statelb) + length(constr.statelb) ...
 +length(constr.initialstatelb) + length(constr.initialstateub) ...
  +length(constr.terminalstatelb) + length(constr.terminalstateub) ~= 6*nx
    error('state constraints are not appropriate size')
elseif length(constr.inputlb)+length(constr.inputub) ~= 2*nu
    error('input constraints are not appropriate size')
end

% Precalculated matrices
[Klqr,P,~] = dlqr(A,B,Q,R); % This P is subject to change
[phi,gamma] = phi_gam_predict(A,B,N);
[omega,psi] = get_omega_psi(Q,P,R,N);
G = 2*(psi+(gamma'*omega*gamma));
F = 2*gamma'*omega*phi;
[W,L,c] = getWLc_own(constr,N,B,gamma,phi);
% [W,L,c] = getWLc(A, B, constr.stateub, constr.statelb, constr.inputub, constr.inputlb, gamma, phi)

% Simulation
xk = [x0 zeros(nx,t)];
uk = zeros(nu,t);
Uk = zeros(nu*N,t);
for k = 1:t
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

figure()
subplot(1,2,1);
stairs(0:t,xk')
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');
xlim([0 t]);
title('Trajectory of states in time')

subplot(1,2,2);
stairs(0:t-1,uk)
xlabel('$k$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xlim([0 t]);
ylim([constr.inputlb-0.1 constr.inputub+0.1]);
title('Inputs in time')
sgtitle('Constrained MPC')

% Comparison with open loop predicted trajecotry
figure()
for k = 1:t
    stairs(k-1:k+N-2,Uk(:,k)) 
    hold on
end
stairs(0:t-1,uk,'k') % MPC control sequence 
xlabel('$k$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xlim([0 t]);
ylim([constr.inputlb-0.1 constr.inputub+0.1]);
title('Open loop predicted trajecotry of input')

% (e) Comparison with saturated unconstrained MPC
% K_mpc = -[eye(nu) zeros(nu,N-1)]*G^-1*F;
% xk = [x0 zeros(nx,t)];
% uk = zeros(nu,t);
% for k = 1:t
%     uk(:,k) = K_mpc*xk(:,k);
%     uk(:,k) = min(max(uk(:,k),constr.inputlb),constr.inputub); % saturation -[1,1]
%     xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
% end
% 
% figure
% subplot(1,2,1);
% stairs(0:t,xk')
% xlabel('$k$','Interpreter','latex');
% ylabel('$x$','Interpreter','latex');
% xlim([0 t]);
% 
% subplot(1,2,2);
% stairs(0:t-1,uk)
% xlabel('$k$','Interpreter','latex');
% ylabel('$u$','Interpreter','latex');
% xlim([0 t]);
% ylim([constr.inputlb-0.1 constr.inputub+0.1]);
% sgtitle('Saturated unconstrained MPC')
