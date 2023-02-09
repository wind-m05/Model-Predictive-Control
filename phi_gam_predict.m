function [phi,gamma] = phi_gam_predict(A,B,N)
phi = zeros(N,1);
gamma = ones(N,N);
gamma = tril(gamma)*B;
for i = 1:N
phi(i) = A^N;
end

for i = 1:N
    for j = 1:N
        gamma(i,j) = A^(i-1) * gamma(i,j);
    end
end


end

