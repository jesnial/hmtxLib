clear variables, close all

%% parameters
alpha = 1e-2;  % singularity shift
N = 6;        % number of points
R0 = 1;
tol = 1e-6;    % tolerance of the ACA method
eta = 3;       % admissibility criterion
Nmin = 3;      % minimum points per cluster

% theta = linspace(-pi/2,pi/2,round(N/3));
% phi = linspace(0,2*pi,round(2*N/3));
% [R,THETA,PHI] = meshgrid(R0,theta,phi);
% P = zeros(numel(R),3);
% [P(:,1),P(:,2),P(:,3)] = sph2cart(R(:),THETA(:),PHI(:));
% P = triSphere(0.2);

% create a grid of N elements in [-1,1] (in 3d for the clustering algorithm)
P = [linspace(-1,1,N).' zeros(N,2)];
% P = [2*rand(N,1)-1 zeros(N,2)];
% P = rand(N,3);

% kernel function: 1/r, protecting the diagonal with alpha
fKern = @(i,j)1./(alpha+pdist2(P(i,:),P(j,:),'euclidean'));

% create the full matrix 
H0 = fKern(1:N,1:N);

% create the H-matrix structure
H = hmtx_cluster(P,'eta',eta,'Nmin',Nmin);

% populate H
H = hmtx_fill(H,tol,fKern);

% plot the H-matrix with the reduced ranks
hmtx_plot(H);
title('H-matrix')

L = hmtx_tril(H);
B = H;

test_lusolve(H,H0,tol);

% X = hmtx_solve(L,B,'left');

%return


% L = hmtx_create('supermatrix',H.irow,H.jcol);
% L = hmtx_copystruct(H);
% L.M{1,2} = hmtx_create('fullmatrix',H.M{1,2}.irow,H.M{1,2}.jcol);
% L.M{1,2}.M = zeros(H.M{1,2}.nrow,H.M{1,2}.ncol);
% L.M{2,1} = hmtx_create('fullmatrix',H.M{2,1}.irow,H.M{2,1}.jcol);
% L.M{2,1}.M = H.M{2,1}.U*H.M{2,1}.V';
% L.M{1,1} = hmtx_create('supermatrix',H.M{1,1}.irow,H.M{1,1}.jcol);
% L.M{1,1}.M{1,1} = hmtx_create('fullmatrix',H.M{1,1}.M{1,1}.irow,H.M{1,1}.M{1,1}.jcol);
% L.M{1,1}.M{1,1}.M = tril(H.M{1,1}.M{1,1}.M);
% L.M{1,1}.M{2,2} = hmtx_create('fullmatrix',H.M{1,1}.M{2,2}.irow,H.M{1,1}.M{2,2}.jcol);
% L.M{1,1}.M{2,2}.M = tril(H.M{1,1}.M{2,2}.M);
% L.M{1,1}.M{2,1} = hmtx_create('fullmatrix',H.M{1,1}.M{2,1}.irow,H.M{1,1}.M{2,1}.jcol);
% L.M{1,1}.M{2,1}.M = H.M{1,1}.M{2,1}.U*H.M{1,1}.M{2,1}.V';
% L.M{1,1}.M{1,2} = hmtx_create('fullmatrix',H.M{1,1}.M{1,2}.irow,H.M{1,1}.M{1,2}.jcol);
% L.M{1,1}.M{1,2}.M = zeros(H.M{1,1}.M{1,2}.nrow,H.M{1,1}.M{1,2}.ncol);
% 
% L.M{2,2} = hmtx_create('supermatrix',H.M{2,2}.irow,H.M{2,2}.jcol);
% L.M{2,2}.M{1,1} = hmtx_create('fullmatrix',H.M{2,2}.M{1,1}.irow,H.M{2,2}.M{1,1}.jcol);
% L.M{2,2}.M{1,1}.M = tril(H.M{2,2}.M{1,1}.M);
% L.M{2,2}.M{2,2} = hmtx_create('fullmatrix',H.M{2,2}.M{2,2}.irow,H.M{2,2}.M{2,2}.jcol);
% L.M{2,2}.M{2,2}.M = tril(H.M{2,2}.M{2,2}.M);
% L.M{2,2}.M{2,1} = hmtx_create('fullmatrix',H.M{2,2}.M{2,1}.irow,H.M{2,2}.M{2,1}.jcol);
% L.M{2,2}.M{2,1}.M = H.M{2,2}.M{2,1}.U*H.M{2,2}.M{2,1}.V';
% L.M{2,2}.M{1,2} = hmtx_create('fullmatrix',H.M{2,2}.M{1,2}.irow,H.M{2,2}.M{1,2}.jcol);
% L.M{2,2}.M{1,2}.M = zeros(H.M{2,2}.M{1,2}.nrow,H.M{2,2}.M{1,2}.ncol);
% 
% Lf = getfield(hmtx_full(L),'M');
% 
% b = ones(H.nrow,1);
% 
% x = hmtx_lsolve(L,b)
% xf = Lf\b
%return
% memory data
[nB,nE] = hmtx_memory(H);
fprintf('Memory compression: %f\n',nB/getfield(whos('H0'),'bytes'));
fprintf('Entry compression: %f\n',nE/numel(H0));

% cluster tree
[qr,qc] = hmtx_getclustertree(H);
a = rand(max(qr),3);
figure, scatter3(P(:,1),P(:,2),P(:,3),80,a(qr,:),'filled');

% create fullmatrix from H-matrix (only for small matrices)
Hf = getfield(hmtx_full(H),'M');

%% check the Frobenius norm of the two matrices
fprintf('Reconstruction accuracy: ')
relnorm = norm(H0-Hf,'fro')/norm(H0,'fro');
if relnorm < tol
    fprintf('\ttest passed: %e < %e\n',relnorm,tol);
else
    fprintf('\t\ttest failed: %e > %e\n',tol);
end

%% transposition
test_transpose(H,H0,tol);

%% scalar multiplication
test_scalar(H,H0,tol);

%% addition and subtraction
test_add(H,H0,tol);
test_sub(H,H0,tol);

%% check matrix multiplication
test_matvect(H,H0,tol);

%% check H-matrix multiplication
test_mult(H,H0,tol);

%% check inverse H-matrix
test_inv(H,H0,tol);

%return
%% test solver and LU decomposition
% create another matrix

% kernel function: 1/r, protecting the diagonal with alpha
beta=3;
gkern = @(i,j)(1./(beta+pdist2(P(i,:),P(j,:), 'euclidean')*5));

Li = gkern(1:N,1:N);
B = hmtx_cluster(P,'eta',eta,'Nmin',Nmin);
B = hmtx_fill(B,tol,gkern);
hmtx_plot(B);
title('B matrix')


test_newsolve(H, H0, B, tol);

test_LU(H, H0,B, tol);


