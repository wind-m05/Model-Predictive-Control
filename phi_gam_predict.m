function [phi,gamma] = phi_gam_predict(A,B,N)
phi = cell(N,1);
gamma = cell(N,N);

zero_vec = zeros(size(A*B));
for j = 1:N
    for i = 1:N
        if i == j
        gamma{i,j} = B;
        elseif i > j
        gamma{i,j} = A^((i-1)-(j-1)) * B;
        else
        gamma{i,j} = zero_vec;
        end
    end
end

for i = 1:N
phi{i} = A^i;
end
phi = cell2mat(phi);
gamma = cell2mat(gamma);

% test1 = eig(phi);
% test2 = eig(gamma);
% if test1<=0 && test2<=0
%     error('Be careful phi and gamma are both negative semi definite')
% end

end

