clear all, close all, clc

%% Simulation parameters
t = 300; 
Ts = 0.1;
N = 20; 
x0 = [10];

%% Constraints
% constr.statelb = [-25;-25;-pi/20;-pi/2;-15]; % v,w,q,theta,h
% constr.stateub =  [25;25;pi/20;pi/2;15]; % v,w,q,theta,h
% constr.initialstatelb = -abs([x0; 0]); % v,w,q,theta,h
% constr.initialstateub = abs([x0; 0]);
% constr.terminalstatelb = [-25;-25;-pi/20;-pi/2;-15]; % Change if needed after monday
% constr.terminalstateub = [25;25;pi/20;pi/2;15];
% constr.inputlb = [-20;-20];
% constr.inputub =  [20;20];
% constr.deltainputlb = [-5;-5];
% constr.deltainputub =  [5;5];

% Model
[A,B,C,D,sys] = modelselect('car1d','discrete',Ts);
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
