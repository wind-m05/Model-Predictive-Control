function [Ccal, Dcal, Ecal, Mcal] = caligraphicMatrices(umin,umax,xmin,xmax,N,n,m)
    Ccal = double.empty;
    for i = 0:(N-1)
        b = [-umin; umax; -xmin; xmax];
        Ccal = [Ccal; b];
    end
    Ccal = [Ccal; -xmin; xmax];
    
    Mcal = double.empty;
    M_i = [zeros(m,n); zeros(m,n); -eye(n); eye(n)];
    for i = 1:(N-1)
        Mcal = blkdiag(Mcal, M_i);
    end
    Mcal = blkdiag(Mcal, [-eye(n); eye(n)]);
    Mcal = [zeros(height(M_i),width(Mcal)); Mcal];
    
    Dcal = [M_i; zeros((height(M_i)*(N-1)+n*2),width(M_i))];
    
    Ecal = double.empty;
    E_i = [-eye(m); eye(m); zeros(n,m); zeros(n,m)];
    for i = 0:(N-1)
        Ecal = blkdiag(Ecal, E_i);
    end
    Ecal = [Ecal; zeros(height([-eye(n); eye(n)]),width(Ecal))];
end