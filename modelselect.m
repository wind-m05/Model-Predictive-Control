function [A,B,C,D] = modelselect(name)
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
   error('Select a model') 
end
end

