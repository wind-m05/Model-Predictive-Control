function [W,L,c,S] = getWLcS(constr,N,B,gamma,phi,T)
%GETWL Summary of this function goes here
%   Detailed explanation goes here
[D,M,E,S,c] = DMESc_constraints(constr,N,B);
W = -D-M*phi;
L = M*gamma+E+S*T;
end