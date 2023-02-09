function [omega,psi] = get_omega_psi(Q,P,R,N)
omega = cell(N,N);
omega_zero_mat = zeros(length(Q),length(Q));
psi_zero_mat = zeros(length(R),length(R));
psi = cell(N,N);
for i = 1:N
    for j = 1:N
        if i == j
            omega{i,j} = Q;
            psi{i,j} = R;
        else
            omega{i,j} = omega_zero_mat;
            psi{i,j} = psi_zero_mat;
        end
    end
end
omega{end,end} = P;
psi = cell2mat(psi);
omega = cell2mat(omega);
end

