clear all, close all, clc
% Simulation parameters
t = 0:1:1500; 
x0 = [3;1];
const.lb = -1; % Lower bound on actuator input 
const.ub = 1; % Upper bound on actuator input
Ts = 0.1;
N = 5;

% Model
[A,B,C,D] = modelselect('car2d');
sys_c = ss(A,B,C,D);
sys_d = c2d(sys_c,Ts);
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;

% Performance, Q square wrt state dim and R square wrt input dim
Q = eye(2);
R = 1;

% Solve LQR and MPC problem 
[Klqr,S,E] = dlqr(Ad,Bd,Q,R);
P = idare(Ad,Bd,Q,R,[],[]);
[phi,gamma] = phi_gam_predict(Ad,Bd,N);
[omega,psi] = get_omega_psi(Q,P,R,N);
G = 2*(psi+(gamma'*omega*gamma));
F = 2*gamma'*omega*phi;
M = zeros(size(B,2),N);
M(1:size(B,2),1:size(B,2)) = eye(size(B,2 ));
Kmpc = M*inv(G)*F; 

% Solve state space iteratively
[y_lqr,x_lqr,t_lqr,u_lqr] = iterativess(sys_d,Klqr,x0,t,const);
[y_mpc,x_mpc,t_mpc,u_mpc] = iterativess(sys_d,Kmpc,x0,t,const);

%% Plotting

figure()
stem(t_lqr(2:end),u_lqr)
hold on
stem(t_mpc(2:end),u_mpc)
title('Actuator input')
legend('lqr','mpc')

figure()
stem(t_lqr,y_lqr(1,:))
hold on
stem(t_lqr,y_lqr(2,:))
stem(t_mpc,y_mpc(1,:))
stem(t_mpc,y_mpc(2,:))
legend('lqr y1','lqr y2','mpc y1','mpc y2')
title('Response')
