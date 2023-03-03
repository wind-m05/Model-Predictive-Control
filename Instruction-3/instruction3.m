clear all, close all, clc

% Simulation parameters
t = 500; 
Ts = 0.1;
N = 10; 
x0 = [3;1];

% Predicted states bounds
constr.statelb = [-10;-10]; 
constr.stateub =  [10;10];
constr.initialstatelb = -abs(x0);
constr.initialstateub = abs(x0);
constr.terminalstatelb = [-0.1;-0.1]; 
constr.terminalstateub = [1;1];
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

% Too computatoinally expensive right now
% %% Feasible sets of inputs and states 
% % Inputs
% Feasibleset_U = Polyhedron('A',L,'b',c+W*x0);
% figure
% if N>3
% plot(projection(Feasibleset_U,[1:3],'vrep')) 
% xlabel('$u_{0|0}$','Interpreter','latex');
% ylabel('$u_{1|0}$','Interpreter','latex');
% zlabel('$u_{2|0}$','Interpreter','latex');
% title('Feasible set of inputs projected onto 3D')
% else
% plot(Feasibleset_U)
% xlabel('$u_{0|0}$','Interpreter','latex');
% ylabel('$u_{1|0}$','Interpreter','latex');
% zlabel('$u_{2|0}$','Interpreter','latex');
% title('Feasible set of inputs')
% end
% % States
% Feasibleset_x0_U = Polyhedron('A',[-W L],'B',c);
% Feasibleset_x0 = projection(Feasibleset_x0_U,1,'vrep');
% % Feasibleset_x0 = projection(Feasibleset_x0_U,1);
% figure
% plot(Feasibleset_x0)
% xlabel('$x_{0|0}$','Interpreter','latex')
% sgtitle('Feasible set of states')

%% Constrained MPC comparing quadprog to activesetsolver

% quadprog
% xk = [x0 zeros(nx,t)];
% uk = zeros(nu,t);
% tk = zeros(1,t);
% for k = 1:t
%     tic;
%     [Uk,fval,exitflag] = quadprog(G,F*xk(:,k),L,c+W*xk(:,k),[],[],[],[],[],[]);
%     tk(k) = toc;
%     if exitflag ~= 1
%         warning('exitflag quadprog = %d\n', exitflag)
%         if exitflag == -2
%             sprintf('Optimization problem is infeasible.')
%         end
%     end
%     uk(:,k) = Uk(1:nu);
%     xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
% end
% 
% figure
% hold on
% xlim([0 t]);
% xlabel('$k$','Interpreter','latex');
% ylabel('$t$','Interpreter','latex');
% sgtitle('Computation time')
% stairs(0:t-1,tk)
% 
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
% 

% (a) MPC Formulation using MPT3
model = LTISystem('A',A,'B',B);
model.x.min = constr.statelb;
model.x.max = constr.stateub;
% or use the following to set state constriant
% model.x.with('setConstraint');
% model.x.setConstraint = Polyhedron('lb', xmin, 'ub', xmax);
model.u.min = constr.inputlb;
model.u.max = constr.inputub;
model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);
mpc = MPCController(model,N);
% to explicit
expmpc = mpc.toExplicit;
% expmpc.partition.plot()

% (b)(c)(d) Comparison with implicit MPC:
% 1. Explicit MPC
t = 300;
xk = [x0 zeros(nx,t)];
uk = zeros(nu,t);
tk = zeros(1,t);
for k = 1:t
    tic;
    uk(:,k) = expmpc.evaluate(xk(:,k));
    tk(k) = toc;
    xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
end

figure
subplot(1,2,1);
stairs(0:t,xk')
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');
xlim([0 t]);

subplot(1,2,2);
stairs(0:t-1,uk)
xlabel('$k$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xlim([0 t]);
sgtitle('Explicit MPC')

f_t = figure;
hold on
xlim([0 t]);
xlabel('$k$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
sgtitle('Computation time')
stairs(0:t-1,tk)

% 2. quadprog
% Compact formulation
xk = [x0 zeros(nx,t)];
uk = zeros(nu,t);
tk = zeros(1,t);
for k = 1:t
    tic;
    [Uk,~,exitflag] = quadprog(G,F*xk(:,k),L,c+W*xk(:,k),[],[],[],[],[],[]);
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
stairs(0:t,xk')
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');
xlim([0 t]);

subplot(1,2,2);
stairs(0:t-1,uk)
xlabel('$k$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
xlim([0 t]);
sgtitle('Online constrained MPC')

figure(f_t);
hold on
xlim([0 t]);
xlabel('$k$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
sgtitle('Computation time')
stairs(0:t-1,tk)
legend('Explicit MPC','Online constrained MPC')
