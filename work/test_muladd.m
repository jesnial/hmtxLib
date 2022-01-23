clear variables, close all

% params
eps = 1e-4;
k12 = 2;
k21 = 3;

% indices
i1 = [1 3 5 7 9]';
i2 = [2 4 6 8]';
j1 = [1 2 3 4 5]';
j2 = [6 7 8 9]';
k1 = [1 2 4 5 7 8]';
k2 = [3 6 9]';

% collect indices;
i = union(i1,i2);
j = union(j1,j2);
k = union(k1,k2);

% matrix A
A = hmtx_create('supermatrix',i,j);
A.eps = eps;
A.M = cell(2,2);

A.M{1,1} = hmtx_create('fullmatrix',i1,j1);
A.M{1,1}.eps = eps;
A.M{1,1}.M = rand(A.M{1,1}.nrow,A.M{1,1}.ncol);

A.M{2,1} = hmtx_create('rkmatrix',i2,j1);
A.M{2,1}.eps = 1e-4;
A.M{2,1}.U = 0*rand(A.M{2,1}.nrow,k21);
A.M{2,1}.V = 0*rand(A.M{2,1}.ncol,k21);
A.M{2,1}.k = k21;

A.M{1,2} = hmtx_create('rkmatrix',i1,j2);
A.M{1,2}.eps = eps;
A.M{1,2}.U = 0*rand(A.M{1,2}.nrow,k12);
A.M{1,2}.V = 0*rand(A.M{1,2}.ncol,k12);
A.M{1,2}.k = k12;

A.M{2,2} = hmtx_create('fullmatrix',i2,j2);
A.M{2,2}.eps = eps;
A.M{2,2}.M = rand(A.M{2,2}.nrow,A.M{2,2}.ncol);

break

% % matrix B
% B = hmtx_create('supermatrix',j,k);
% B.eps = eps;
% B.M = cell(2,2);
% 
% B.M{1,1} = hmtx_create('fullmatrix',j1,k1);
% B.M{1,1}.eps = eps;
% B.M{1,1}.M = rand(B.M{1,1}.nrow,B.M{1,1}.ncol);
% 
% B.M{2,1} = hmtx_create('rkmatrix',j2,k1);
% B.M{2,1}.eps = eps;
% B.M{2,1}.U = rand(B.M{2,1}.nrow,k21);
% B.M{2,1}.V = rand(B.M{2,1}.ncol,k21);
% B.M{2,1}.k = k21;
% 
% B.M{1,2} = hmtx_create('rkmatrix',j1,k2);
% B.M{1,2}.eps = eps;
% B.M{1,2}.U = rand(B.M{1,2}.nrow,k12);
% B.M{1,2}.V = rand(B.M{1,2}.ncol,k12);
% B.M{1,2}.k = k12;
% 
% B.M{2,2} = hmtx_create('fullmatrix',j2,k2);
% B.M{2,2}.eps = eps;
% B.M{2,2}.M = rand(B.M{2,2}.nrow,B.M{2,2}.ncol);

% % matrix B
% B = hmtx_create('supermatrix',j,k);
% B.eps = eps;
% B.M = cell(2,2);
% 
% B.M{1,1} = hmtx_create('fullmatrix',j1,k1);
% B.M{1,1}.eps = eps;
% B.M{1,1}.M = rand(B.M{1,1}.nrow,B.M{1,1}.ncol);
% 
% B.M{2,1} = hmtx_create('fullmatrix',j2,k1);
% B.M{2,1}.eps = eps;
% B.M{2,1}.M = rand(B.M{2,1}.nrow,B.M{2,1}.ncol);
% 
% B.M{1,2} = hmtx_create('rkmatrix',j1,k2);
% B.M{1,2}.eps = eps;
% B.M{1,2}.U = rand(B.M{1,2}.nrow,k12);
% B.M{1,2}.V = rand(B.M{1,2}.ncol,k12);
% B.M{1,2}.k = k12;
% 
% B.M{2,2} = hmtx_create('fullmatrix',j2,k2);
% B.M{2,2}.eps = eps;
% B.M{2,2}.M = rand(B.M{2,2}.nrow,B.M{2,2}.ncol);

% matrix C
C = hmtx_create('supermatrix',i,k);
C.eps = eps;
C.M = cell(2,2);

C.M{1,1} = hmtx_create('fullmatrix',i1,k1);
C.M{1,1}.eps = eps;
C.M{1,1}.M = rand(C.M{1,1}.nrow,C.M{1,1}.ncol);

C.M{2,1} = hmtx_create('rkmatrix',i2,k1);
C.M{2,1}.eps = 1e-4;
C.M{2,1}.U = rand(C.M{2,1}.nrow,k21);
C.M{2,1}.V = rand(C.M{2,1}.ncol,k21);
C.M{2,1}.k = k21;

C.M{1,2} = hmtx_create('rkmatrix',i1,k2);
C.M{1,2}.eps = eps;
C.M{1,2}.U = rand(C.M{1,2}.nrow,k12);
C.M{1,2}.V = rand(C.M{1,2}.ncol,k12);
C.M{1,2}.k = k12;

C.M{2,2} = hmtx_create('fullmatrix',i2,k2);
C.M{2,2}.eps = eps;
C.M{2,2}.M = rand(C.M{2,2}.nrow,C.M{2,2}.ncol);

hmtx_plot(A);
hmtx_plot(B);
hmtx_plot(C);

Af = hmtx_full(A);
Bf = hmtx_full(B);
Cf = hmtx_full(C);

C1 = Cf.M+Af.M*Bf.M;

C = hmtx_muladd2(C,A,B);
Cf = hmtx_full(C);
C2 = Cf.M;

norm(C1-C2,'fro')

