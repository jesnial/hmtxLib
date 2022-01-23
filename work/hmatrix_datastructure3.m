clear variables, close all
M = 2^4;
N = 2^4;
% vector
x = rand(M,1);

% matrix
% A = cell(2,2);
A.type = 'supermatrix';
A.m = int32(1:M)';
A.n = int32(1:N)';
A.M = cell(2,2);

A.M{1,1}.type = 'fullmatrix';
A.M{1,1}.m = int32(1:M/2)';
A.M{1,1}.n = int32(1:N/2)';
A.M{1,1}.M = rand(A.M{1,1}.m,A.M{1,1}.n);

A.M{1,2}.type = 'fullmatrix';
A.M{1,2}.m = M/2;
A.M{1,2}.n = N/2;
A.M{1,2}.M = rand(A.M{1,2}.m,A.M{1,2}.n);

A.M{2,1}.type = 'rkmatrix';
A.M{2,1}.m = M/2;
A.M{2,1}.n = N/2;
A.M{2,1}.r = 3;
A.M{2,1}.U = rand(A.M{2,1}.m,A.M{2,1}.r);
A.M{2,1}.VT = rand(A.M{2,1}.r,A.M{2,1}.n);

A.M{2,2}.type = 'supermatrix';
A.M{2,2}.m = M/2;
A.M{2,2}.n = N/2;
A.M{2,2}.M = cell(2,2);

A.M{2,2}.M{1,1}.type = 'fullmatrix';
A.M{2,2}.M{1,1}.m = M/4;
A.M{2,2}.M{1,1}.n = N/4;
A.M{2,2}.M{1,1}.M = rand(A.M{2,2}.M{1,1}.m,A.M{2,2}.M{1,1}.n);
A.M{2,2}.M{1,2}.type = 'fullmatrix';
A.M{2,2}.M{1,2}.m = M/4;
A.M{2,2}.M{1,2}.n = N/4;
A.M{2,2}.M{1,2}.M = rand(A.M{2,2}.M{1,1}.m,A.M{2,2}.M{1,1}.n);
A.M{2,2}.M{2,1}.type = 'fullmatrix';
A.M{2,2}.M{2,1}.m = M/4;
A.M{2,2}.M{2,1}.n = N/4;
A.M{2,2}.M{2,1}.M = rand(A.M{2,2}.M{1,1}.m,A.M{2,2}.M{1,1}.n);

A.M{2,2}.M{2,2}.type = 'supermatrix';
A.M{2,2}.M{2,2}.m = M/4;
A.M{2,2}.M{2,2}.n = N/4;
A.M{2,2}.M{2,2}.M = cell(2,2);
A.M{2,2}.M{2,2}.M{1,1}.type = 'fullmatrix';
A.M{2,2}.M{2,2}.M{1,1}.m = M/8;
A.M{2,2}.M{2,2}.M{1,1}.n = N/8;
A.M{2,2}.M{2,2}.M{1,1}.M = rand(A.M{2,2}.M{2,2}.M{1,1}.m,A.M{2,2}.M{2,2}.M{1,1}.n);
A.M{2,2}.M{2,2}.M{1,2}.type = 'fullmatrix';
A.M{2,2}.M{2,2}.M{1,2}.m = M/8;
A.M{2,2}.M{2,2}.M{1,2}.n = N/8;
A.M{2,2}.M{2,2}.M{1,2}.M = rand(A.M{2,2}.M{2,2}.M{1,1}.m,A.M{2,2}.M{2,2}.M{1,1}.n);
A.M{2,2}.M{2,2}.M{2,1}.type = 'fullmatrix';
A.M{2,2}.M{2,2}.M{2,1}.m = M/8;
A.M{2,2}.M{2,2}.M{2,1}.n = N/8;
A.M{2,2}.M{2,2}.M{2,1}.M = rand(A.M{2,2}.M{2,2}.M{1,1}.m,A.M{2,2}.M{2,2}.M{1,1}.n);
A.M{2,2}.M{2,2}.M{2,2}.type = 'fullmatrix';
A.M{2,2}.M{2,2}.M{2,2}.m = M/8;
A.M{2,2}.M{2,2}.M{2,2}.n = N/8;
A.M{2,2}.M{2,2}.M{2,2}.M = rand(A.M{2,2}.M{2,2}.M{1,1}.m,A.M{2,2}.M{2,2}.M{1,1}.n);

% assemble full matrix
Am = fullmatrix(A);

tic
b = Am*x;
toc

tic
b2 = hmtx_Ax(A,x);
toc

h = plothmatrix(A);
