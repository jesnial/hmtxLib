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




B = hmtx_create('supermatrix',1:n1+n2,1:m1+m2);
B.eps = 1e-4;
B.M = cell(2,2);

B.M{1,1} = hmtx_create('fullmatrix',1:2:n1+n2,1:2:m1+m2);
B.M{1,1}.type = 'fullmatrix';
B.M{1,1}.eps = 1e-4;
B.M{1,1}.M = rand(B.M{1,1}.nrow,B.M{1,1}.ncol)+eye(B.M{1,1}.nrow,B.M{1,1}.ncol);

B.M{2,2} = hmtx_create('fullmatrix',2:2:n1+n2,2:2:m1+m2);
B.M{2,2}.eps = 1e-4;
B.M{2,2}.M = rand(B.M{2,2}.nrow,B.M{2,2}.ncol)+eye(B.M{2,2}.nrow,B.M{2,2}.ncol);

B.M{2,1} = hmtx_create('rkmatrix',2:2:n1+n2,1:2:m1+m2);
B.M{2,1}.eps = 1e-4;
B.M{2,1}.U = rand(B.M{2,1}.nrow,k12);
B.M{2,1}.V = rand(B.M{2,1}.ncol,k12);
B.M{2,1}.k = k12;


B.M{1,2} = hmtx_create('rkmatrix',1:2:n1+n2,2:2:m1+m2);
B.M{1,2}.eps = 1e-4;
B.M{1,2}.U = rand(B.M{1,2}.nrow,k21);
B.M{1,2}.V = rand(B.M{1,2}.ncol,k21);
B.M{1,2}.k = k21;


hmtx_plot(A);
Af = hmtx_full(A);

B = hmtx_add(A,A,1);
