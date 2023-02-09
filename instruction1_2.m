clear all, close all, clc
% First order model
A = 0;
B = 1;
C = 1;
Ts = 0.1;
sys_c = ss(A,B,C,0);
sys_d = c2d(sys_c,Ts);
Q = 1;
R = 100; % R is 1 saturates the input to 1, R = 100 saturates the input to 0.1, takes 10 times as long
[K,S,e] = dlqr(sys_d.A,sys_d.B,Q,R);
sys_cl = ss(sys_d.A-sys_d.B*K,sys_d.B,sys_d.C,0,Ts);
[y,t,x] = initial(sys_cl,1);
figure()
stem(t,y)
figure()
stem(t,-K*x)