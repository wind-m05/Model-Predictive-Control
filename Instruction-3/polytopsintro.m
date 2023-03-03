%% Polytops intro
% 2D polytop
A = [-1,0;0.5,0;0,-0.5;0,0.25;1,1];
b = [1;1;1;1;1];
P = Polyhedron(A,b);
P.plot()
% 3D polytop and projection
A = [1,0,0;
     -1,0,0;
      0,1,0;
      0,-1,0;
      0,0,1;
      0,0,-1];
b=[4;4;1;1;2;2];
P2 = Polyhedron(A,b);
figure
P2.plot()
P12 = projection(P2,[1])
figure
P12.plot()
%% Example mpt3
% p1 = Polyhedron(L,c+W*x0)
A = [1,1;0,1];
B = [1;0.5];
model = LTISystem('A',A,'B',B);
model.x.min = [-5;-5];
model.x.max = [5;5];
model.u.min = -1;
model.u.max = 1;
Q = eye(2);
R = 1;
model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);
N = 5;
mpc = MPCController(model,N);
x = [3;-1];
% Explicit solution
expmpc = mpc.toExplicit();
u = expmpc.evaluate(x);
expmpc.partition.plot()