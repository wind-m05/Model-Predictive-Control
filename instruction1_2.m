clear all, close all, clc
% First order model
A = [0 1;
     0 0];
B = [0;1];
C = eye(2);
Ts = 0.1;
x0 = [3;1];
N = 5;
sys_c = ss(A,B,C,0);
sys_d = c2d(sys_c,Ts);
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Q = eye(2);
R = 1; % R is 1 saturates the input to 1, R = 100 saturates the input to 0.1, takes 10 times as long
[K,S,E] = dlqr(Ad,Bd,Q,R);
sys_cl = ss(Ad-Bd*K,Bd,Cd,0,Ts);
[y,t,x] = initial(sys_cl,x0);
P = idare(Ad,Bd,Q,R,[],[]);
[phi,gamma] = phi_gam_predict(Ad,Bd,N);
[omega,psi] = get_omega_psi(Q,P,R,N);
G = 2*(psi+gamma'*omega*gamma);
F = 2*gamma'*omega*psi;
Binpre = eye(length(B));
M = zeros(size(G));
M(1:length(B),1:length(B)) = eye(length(B));
Kmpc = -M*(inv(G)*F);

% sys_cl = ss(Ad-Bd*Kmpc,Bd,Cd,0,Ts);
% [y,t,x] = initial(sys_cl,x0);
% figure()
% stem(t,y)
% figure()
% stem(t,-K*x)


%% MPC contraints
% t = 1:0.1:10;
% xold = 1;
% 
% for i = 1:length(t)
% umpc = Kmpc*xold;
% xnew = A*xold+B*umpc;
% ynew{i} = C*xnew;
% xold = xnew;
% end
% ynew = cell2mat(ynew);
% 
% figure()
% plot(t,ynew)
