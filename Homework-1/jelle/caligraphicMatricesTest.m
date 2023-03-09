function [CcalTest, DcalTest, EcalTest, McalTest] = caligraphicMatricesTest(umin,umax,xmin,xmax,N,n,m)
    CcalTest = zeros((height(umin)+height(umax)+height(xmin)+height(xmax))*N+height(xmin)+height(xmax),1);
    DcalTest = zeros((height(umin)+height(umax)+height(xmin)+height(xmax))*N+height(xmin)+height(xmax),n);
    McalTest = zeros((height(umin)+height(umax)+height(xmin)+height(xmax))*N+height(xmin)+height(xmax),N*n);
    EcalTest = zeros((height(umin)+height(umax)+height(xmin)+height(xmax))*N+height(xmin)+height(xmax),N*m);
end