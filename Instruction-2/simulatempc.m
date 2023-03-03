function [y,x,t,uhist] = simulatempc(sys,x0,t,const)
x{1} = x0;
y{1} = sys.C*x{1};

for i = 1:length(t)-1
    [U,~,exitflag] = quadprog(G,F*xk(:,k),L,c+W*xk(:,k),[],[],[],[],[],opt);
    if exitflag ~= 1
        warning('exitflag quadprog = %d\n', exitflag)
        if exitflag == -2
            sprintf('Optimization problem is infeasible.')
        end
    end
    u = constrain(-K*y{i},const.lb,const.ub);
    uhist{i} = u;
    x{i+1} = sys.A*x{i}+sys.B*u;
    y{i+1} = sys.C*x{i+1};
end
y = cell2mat(y);
x = cell2mat(x);
uhist = cell2mat(uhist);
end


