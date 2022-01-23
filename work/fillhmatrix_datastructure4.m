clear variables, close all
M = 2^4;
N = 2^4;
% vector
x = rand(M,1);

B = [rand(M/2,3); 2+rand(M/2,3)];

% matrix
% A = cell(2,2);
A.type = 'supermatrix';
A.m = M;
A.n = N;
A.im = int32(1:M)';
A.in = int32(1:N)';
A.M = cell(2,2);

A.M{1,1}.type = 'fullmatrix';
A.M{1,1}.m = M/2;
A.M{1,1}.n = N/2;
A.M{1,1}.im = int32(1:M/2)';
A.M{1,1}.in = int32(1:N/2)';
%A.M{1,1}.M = rand(A.M{1,1}.m,A.M{1,1}.n);

A.M{1,2}.type = 'fullmatrix';
A.M{1,2}.m = M/2;
A.M{1,2}.n = N/2;
A.M{1,2}.im = int32(1:M/2)';
A.M{1,2}.in = int32(N/2+1:N)';
%A.M{1,2}.M = rand(A.M{1,2}.m,A.M{1,2}.n);

A.M{2,1}.type = 'rkmatrix';
A.M{2,1}.m = M/2;
A.M{2,1}.n = N/2;
A.M{2,1}.im = int32(M/2+1:N)';
A.M{2,1}.in = int32(1:N/2)';
%A.M{2,1}.r = 3;
%A.M{2,1}.U = rand(A.M{2,1}.m,A.M{2,1}.r);
%A.M{2,1}.VT = rand(A.M{2,1}.r,A.M{2,1}.n);

A.M{2,2}.type = 'fullmatrix';
A.M{2,2}.m = M/2;
A.M{2,2}.n = N/2;
A.M{2,2}.im = int32(M/2+1:N)';
A.M{2,2}.in = int32(N/2+1:N)';

% A.M{2,2}.type = 'supermatrix';
% A.M{2,2}.m = M/2;
% A.M{2,2}.n = N/2;
% A.M{2,2}.im = int32(M/2+1:N)';
% A.M{2,2}.in = int32(N/2+1:N)';
% A.M{2,2}.M = cell(2,2);
% 
% A.M{2,2}.M{1,1}.type = 'fullmatrix';
% A.M{2,2}.M{1,1}.m = M/4;
% A.M{2,2}.M{1,1}.n = N/4;
% A.M{2,2}.M{1,1}.im = int32(M/2+1:3*M/4);
% A.M{2,2}.M{1,1}.in = int32(N/2+1:3*N/4);
% %A.M{2,2}.M{1,1}.M = rand(A.M{2,2}.M{1,1}.m,A.M{2,2}.M{1,1}.n);
% A.M{2,2}.M{1,2}.type = 'fullmatrix';
% A.M{2,2}.M{1,2}.m = M/4;
% A.M{2,2}.M{1,2}.n = N/4;
% A.M{2,2}.M{1,2}.im = int32(M/2+1:3*M/4);
% A.M{2,2}.M{1,2}.in = int32(3*N/4+1:N);
% %A.M{2,2}.M{1,2}.M = rand(A.M{2,2}.M{1,1}.m,A.M{2,2}.M{1,1}.n);
% A.M{2,2}.M{2,1}.type = 'fullmatrix';
% A.M{2,2}.M{2,1}.m = M/4;
% A.M{2,2}.M{2,1}.n = N/4;
% A.M{2,2}.M{2,1}.im = int32(M/2+1:3*M/4);
% A.M{2,2}.M{2,1}.in = int32(3*N/4+1:N);
% %A.M{2,2}.M{2,1}.M = rand(A.M{2,2}.M{1,1}.m,A.M{2,2}.M{1,1}.n);
% 
% A.M{2,2}.M{2,2}.type = 'supermatrix';
% A.M{2,2}.M{2,2}.m = M/4;
% A.M{2,2}.M{2,2}.n = N/4;
% A.M{2,2}.M{2,2}.im = int32(3*M/4+1:M);
% A.M{2,2}.M{2,2}.in = int32(3*N/4+1:N);
% A.M{2,2}.M{2,2}.M = cell(2,2);
% A.M{2,2}.M{2,2}.M{1,1}.type = 'fullmatrix';
% A.M{2,2}.M{2,2}.M{1,1}.m = M/8;
% A.M{2,2}.M{2,2}.M{1,1}.n = N/8;
% A.M{2,2}.M{2,2}.M{1,1}.im = (3*M/4+1:7*M/8);
% A.M{2,2}.M{2,2}.M{1,1}.in = (3*N/4+1:7*N/8);
% %A.M{2,2}.M{2,2}.M{1,1}.M = rand(A.M{2,2}.M{2,2}.M{1,1}.m,A.M{2,2}.M{2,2}.M{1,1}.n);
% A.M{2,2}.M{2,2}.M{1,2}.type = 'fullmatrix';
% A.M{2,2}.M{2,2}.M{1,2}.m = M/8;
% A.M{2,2}.M{2,2}.M{1,2}.n = N/8;
% A.M{2,2}.M{2,2}.M{1,2}.im = (3*M/4+1:7*M/8);
% A.M{2,2}.M{2,2}.M{1,2}.in = (7*N/8+1:N);
% %A.M{2,2}.M{2,2}.M{1,2}.M = rand(A.M{2,2}.M{2,2}.M{1,1}.m,A.M{2,2}.M{2,2}.M{1,1}.n);
% A.M{2,2}.M{2,2}.M{2,1}.type = 'fullmatrix';
% A.M{2,2}.M{2,2}.M{2,1}.m = M/8;
% A.M{2,2}.M{2,2}.M{2,1}.n = N/8;
% A.M{2,2}.M{2,2}.M{2,1}.im = (7*M/8+1:M);
% A.M{2,2}.M{2,2}.M{2,1}.in = (3*N/4+1:7*N/8);
% %A.M{2,2}.M{2,2}.M{2,1}.M = rand(A.M{2,2}.M{2,2}.M{1,1}.m,A.M{2,2}.M{2,2}.M{1,1}.n);
% A.M{2,2}.M{2,2}.M{2,2}.type = 'fullmatrix';
% A.M{2,2}.M{2,2}.M{2,2}.m = M/8;
% A.M{2,2}.M{2,2}.M{2,2}.n = N/8;
% A.M{2,2}.M{2,2}.M{2,2}.im = (7*M/8+1:M);
% A.M{2,2}.M{2,2}.M{2,2}.in = (7*N/8+1:N);
% %A.M{2,2}.M{2,2}.M{2,2}.M = rand(A.M{2,2}.M{2,2}.M{1,1}.m,A.M{2,2}.M{2,2}.M{1,1}.n);

fkern = @(irow,jcol)one_dividedby_r(irow,jcol,B,B);

eps = 1e-6;
tic
A = fillhmatrix(A,fkern,eps);
toc
tic
Af = fkern(1:M,1:N);
toc


% assemble full matrix
Am = fullmatrix(A);

tic
b = Am*x;
toc

tic
b2 = hmtx_Ax(A,x);
toc

h = plothmatrix(A);
