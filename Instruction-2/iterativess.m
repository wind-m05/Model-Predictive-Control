function [y,x,t,uhist] = iterativess(sys,K,x0,t,const)
x{1} = x0;
y{1} = sys.C*x{1};

for i = 1:length(t)-1
    u = constrain(-K*y{i},const.lb,const.ub);
    uhist{i} = u;
    x{i+1} = sys.A*x{i}+sys.B*u;
    y{i+1} = sys.C*x{i+1};
end
y = cell2mat(y);
x = cell2mat(x);
uhist = cell2mat(uhist);
end

