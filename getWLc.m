function [W,L,c] = getWLc(constr,N,Bd,gamma,phi)
%GETWL Summary of this function goes here
%   Detailed explanation goes here
[D,M,E,c] = DMEc_constraints(constr,N,Bd);
W = -D-M*phi;
L = M*gamma+E;
end

