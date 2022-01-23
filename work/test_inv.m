clear variables, close all
n1 = 5;
n2 = 6;
m1 = 5;
m2 = 6;
k12 = 3;
k21 = 2;

A = hmtx_create('supermatrix',1:n1+n2,1:m1+m2);
A.eps = 1e-4;
A.M = cell(2,2);

A.M{1,1} = hmtx_create('fullmatrix',1:2:n1+n2,1:2:m1+m2);
A.M{1,1}.type = 'fullmatrix';
A.M{1,1}.eps = 1e-4;
A.M{1,1}.M = rand(A.M{1,1}.nrow,A.M{1,1}.ncol)+eye(A.M{1,1}.nrow,A.M{1,1}.ncol);

A.M{2,2} = hmtx_create('fullmatrix',2:2:n1+n2,2:2:m1+m2);
A.M{2,2}.eps = 1e-4;
A.M{2,2}.M = rand(A.M{2,2}.nrow,A.M{2,2}.ncol)+eye(A.M{2,2}.nrow,A.M{2,2}.ncol);

A.M{1,2} = hmtx_create('rkmatrix',1:2:n1+n2,2:2:m1+m2);
A.M{1,2}.eps = 1e-4;
A.M{1,2}.U = rand(A.M{1,2}.nrow,k12);
A.M{1,2}.V = rand(A.M{1,2}.ncol,k12);
A.M{1,2}.k = k12;


A.M{2,1} = hmtx_create('rkmatrix',2:2:n1+n2,1:2:m1+m2);
A.M{2,1}.eps = 1e-4;
A.M{2,1}.U = rand(A.M{2,1}.nrow,k21);
A.M{2,1}.V = rand(A.M{2,1}.ncol,k21);
A.M{2,1}.k = k21;

% hmtx_plot(A);
% Af = hmtx_full(A);
% Ainv = hmtx_inv(A);
% Ainvf = hmtx_full(Ainv);
% Afinv = inv(Af.M);
% norm(Ainvf.M)
% norm(Afinv)
% norm(Ainvf.M-Afinv)

B = hmtx_create('rkmatrix',A.irow,A.jcol);
B.eps = 1e-6;
B.k = 2;
B.U = rand(A.nrow,B.k);
B.V = rand(A.ncol,B.k);

% full matrices
Af = hmtx_full(A);
Bf = hmtx_full(B);
C1 = Af.M+Bf.M;
C = hmtx_add(Af,Bf,'+');
Cf = hmtx_full(C);
C2 = Cf.M;

B2 = rk2super(B,A.M{1,1}.irow,A.M{1,1}.jcol);
C = hmtx_add(A,B2,'+');

Cff = hmtx_full(C);
C3 = Cff.M;

