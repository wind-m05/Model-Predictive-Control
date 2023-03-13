function [phi, gamma] = predictionModel(A,B,N,n,m)
    phi = double.empty;
    gamma = zeros(N*n,N*m);
    for i = 1:N
        phi = [phi; A^i];
        for j = 1:i
            gamma((i-1)*n+1:n*i,(j-1)*m+1:j*m) = A^(i-j)*B;
        end
    end
end