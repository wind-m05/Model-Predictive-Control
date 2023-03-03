function [A,B,C,D,sys] = modelselect(name,cod,Ts)
if strcmp(name,'car1d')
    A = 0;
    B = 1;
    C = 1;
    D = 0;
elseif strcmp(name,'car2d')
    A = [0 1;
        0 0];
    B = [0;-1];
    C = eye(2);
    D = 0;
else
   error('Select a valid model') 
end

sys_temp = ss(A,B,C,D);

if strcmp(cod,'discrete')
sys = c2d(sys_temp,Ts,'zoh');
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;
elseif strcmp(cod,'continuous')
A = sys_temp.A;
B = sys_temp.B;
C = sys_temp.C;
D = sys_temp.D;
sys = sys_temp;
end
end

