function T = T_delta(N,nu)
%T_DELTA Summary of this function goes here
%   Detailed explanation goes here
T = cell(N,N)
for i = 1:N
    for j = 1:N
        if i == j
            T{i,j} = eye(nu);
        elseif i == j+1
            T{i,j} = -eye(nu);
        else
            T{i,j} = zeros(nu);
        end
    end
end
T = cell2mat(T);
end

