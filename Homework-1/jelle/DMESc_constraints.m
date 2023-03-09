function [D,M,E,S,c] = DMESc_constraints(constr,N,B_sys)
%DME_CONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here
nx = size(B_sys,1); 
nu = size(B_sys,2);
ncx = 2*length(constr.statelb);
ncu = 2*(length(constr.inputlb)+length(constr.deltainputlb));

m{1} = zeros(ncu,nx);    
m{2} = [-eye(nx);0 1 0 -7.74];
m{3} = [eye(nx);0 -1 0 7.74];

e{1} = -eye(nu);  
e{2} = eye(nu);  
e{3} = zeros(nx,nu);  
e{4} = zeros(ncx,nu); 

s{1} = zeros(nx,nu);
s{2} = -eye(nu);
s{3} = eye(nu);
s{4} = zeros(ncx,nu);

% bi
b{1} = -constr.inputlb;
b{2} = constr.inputub;
b{3} = -constr.deltainputlb;
b{4} = constr.deltainputub;
b{5} = -constr.statelb;
b{6} = constr.stateub;
% Define other b's for b0 and bN ?? 

m_mat = cell2mat(m');
e_mat = cell2mat(e');
s_mat = cell2mat(s');
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

% Curly S
for i = 1:N+1
    for j = 1:N 
        if i == j
        S{i,j} = s_mat;
        else
        S{i,j} = zeros(size(s_mat));
        end
    end
end

% c (change if you also need specific b0 and bN)
for i = 1:N+1
B{i,1} = b_mat;
end

D = cell2mat(D);
M = cell2mat(M);
E = cell2mat(E);
S = cell2mat(S);
c = cell2mat(B);
end