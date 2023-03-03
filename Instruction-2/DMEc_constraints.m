function [D,M,E,c] = DMEc_constraints(constr,N,Bd)
%DME_CONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here
nx = size(Bd,1);
nu = size(Bd,2);
% m = cell()
m{1} = zeros(nu,nx);   
m{2} = zeros(nu,nx);  
m{3} = -eye(nx);
m{4} = eye(nx);
e{1} = -eye(nu);  
e{2} = eye(nu);  
e{3} = zeros(nx,nu);  
e{4} = zeros(nx,nu);  
b{1} = -constr.inputlb;
b{2} = constr.inputub;
b{3} = -constr.statelb;
b{4} = constr.stateub;
m_mat = cell2mat(m');
e_mat = cell2mat(e');
b_mat = cell2mat(b');
% Curly D
D{1} = m_mat;
for i = 2:N+1
    D{i,1} = zeros(size(m_mat));
end

% Curly M
for i = 1:N+1
    for j = 1:N   
        if i == j +1
        M{i,j} = m_mat;
        else
        M{i,j} = zeros(size(m_mat)); 
        end
    end
end

% Curly E
for i = 1:N+1
    for j = 1:N 
        if i == j
        E{i,j} = e_mat;
        else
        E{i,j} = zeros(size(e_mat));
        end
    end
end

for i = 1:N+1
B{i,1} = b_mat;
end

D = cell2mat(D);
M = cell2mat(M);
E = cell2mat(E);
c = cell2mat(B);
end
