clear all, close all, clc
%% 2.1 a)
Ts = 1;
% Model
[A,B,C,D] = modelselect('car1d');
sys_c = ss(A,B,C,D);
sys_d = c2d(sys_c,Ts);
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;

%% Optimization
options = optimoptions('quadprog','Display','iter');
H = 2;
f = 5;
A = [];
b = [];
Aeq = [];
beq = [];
lb = -1;
ub = 1;
x0 = [];
x = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,[])

%% 2.1 b)
Ts = 1;

% Model
[A,B,C,D] = modelselect('car2d');
sys_c = ss(A,B,C,D);
sys_d = c2d(sys_c,Ts);
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;

%% Optimization
options = optimoptions('quadprog','Display','iter');
H = 2.25;
f = 9;
A = [];
b = [];
Aeq = [];
beq = [];
lb = -1;
ub = 1;
x0 = [];
x = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,[])