clear all, close all, clc

% Predicted states bounds
constr.statelb = [-1;-1]; 
constr.stateub =  [1;1];
% Current state bound 
constr.initialstatelb = -0.5;
constr.initialstateub = 0.5;
% Terminal state bound
constr.terminalstatelb = -3;
constr.terminalstateub = 3;
% Input bounds
constr.inputlb = -2;
constr.inputub =  2;

% Simulation parameters
t = 0:1:300; 
x0 = [3;1];
N = 5; 
Ts = 0.1;

% Model
[A,B,C,D] = modelselect('car2d');
sys_c = ss(A,B,C,D);
sys_d = c2d(sys_c,Ts,'zoh');
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Q = eye(2);
R = 1;
P = idare(Ad,Bd,Q,R,[],[]);
[phi,gamma] = phi_gam_predict(Ad,Bd,N);
[omega,psi] = get_omega_psi(Q,P,R,N);
G = 2*(psi+(gamma'*omega*gamma));
F = 2*gamma'*omega*phi;
[W,L,c] = getWLc(constr,N,Bd,gamma,phi);
opt =  optimoptions('quadprog','Display','off')
[y,x,t,u] = simulatempc(sys_d,Klqr,x0,t,const)
xk = [x0 zeros(nx,k_sim)]
% Simulation
% % (c) Closed-loop trajectory
% k_sim = 500;
% xk = [x0 zeros(nx,k_sim)];
% uk = zeros(nu,k_sim);
% Uk = zeros(nu*N,k_sim);
% opt =  optimoptions('quadprog','Display','off');
% warning('off','optim:quadprog:HessianNotSym')
% for k = 1:k_sim
%     [U,~,exitflag] = quadprog(G,F*xk(:,k),L,c+W*xk(:,k),[],[],[],[],[],opt);
%     if exitflag ~= 1
%         warning('exitflag quadprog = %d\n', exitflag)
%         if exitflag == -2
%             sprintf('Optimization problem is infeasible.')
%         end
%     end
%     Uk(:,k) = U;
%     uk(:,k) = U(1:nu);
%     xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
% end
% 
% figure
% subplot(1,2,1);
% stairs(0:k_sim,xk')
% xlabel('$k$','Interpreter','latex');
% ylabel('$x$','Interpreter','latex');
% xlim([0 k_sim]);
% 
% subplot(1,2,2);
% stairs(0:k_sim-1,uk)
% xlabel('$k$','Interpreter','latex');
% ylabel('$u$','Interpreter','latex');
% xlim([0 k_sim]);
% ylim([umin-0.1 umax+0.1]);
% sgtitle('Constrained MPC')
% 
% % (d) Comparison with open loop predicted trajecotry
% % Follow the previous section
% figure
% hold on
% for k = 1:k_sim
%     stairs(k-1:k+N-2,Uk(:,k))
% end
% stairs(0:k_sim-1,uk,'k')
% xlabel('$k$','Interpreter','latex');
% ylabel('$u$','Interpreter','latex');
% xlim([0 k_sim]);
% ylim([umin-0.1 umax+0.1]);
% title('Open loop predicted trajecotry of input')
% 
% % (e) Comparison with saturated unconstrained MPC
% K_mpc = -[eye(nu) zeros(nu,N-1)]*G^-1*F;
% k_sim = 500;
% xk = [x0 zeros(nx,k_sim)];
% uk = zeros(nu,k_sim);
% for k = 1:k_sim
%     uk(:,k) = K_mpc*xk(:,k);
%     uk(:,k) = min(max(uk(:,k),umin),umax); % saturation -[1,1]
%     xk(:,k+1) = A*xk(:,k)+B*uk(:,k);
% end
% 
% figure
% subplot(1,2,1);
% stairs(0:k_sim,xk')
% xlabel('$k$','Interpreter','latex');
% ylabel('$x$','Interpreter','latex');
% xlim([0 k_sim]);
% 
% subplot(1,2,2);
% stairs(0:k_sim-1,uk)
% xlabel('$k$','Interpreter','latex');
% ylabel('$u$','Interpreter','latex');
% xlim([0 k_sim]);
% ylim([umin-0.1 umax+0.1]);
% sgtitle('Saturated unconstrained MPC')