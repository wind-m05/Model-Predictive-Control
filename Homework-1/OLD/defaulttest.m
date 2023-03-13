function [x0,constr,Q,R] = defaulttest(nx,nu)
%DEFAULTTEST Summary of this function goes here
%   Detailed explanation goes here
x0 = ones(nx,1);
constr.statelb = -ones(nx,1);
constr.stateub =  ones(nx,1);
constr.initialstatelb = -ones(nx,1);
constr.initialstateub = ones(nx,1);
constr.terminalstatelb = -ones(nx,1);
constr.terminalstateub = ones(nx,1);
constr.inputlb = -ones(nu,1);
constr.inputub =  ones(nu,1);
Q = 1*eye(nx);
R = 1*eye(nu);
end

