clear variables, close all

% matrix generation
alpha = 1e-2;
N = 10;
tol = 1e-8;  % tolerance of the ACA method
eta = 3;
Nmin = 5;

% create a grid of N elements in [-1,1]. The clustering algorithm works
% with 3d points
P = [linspace(-1,1,N).' zeros(N,2)];
% P = [2*rand(N,1)-1 zeros(N,2)];
% P = rand(N,3);

% kernel function: 1/r, protecting the diagonal with alpha
% fkern = @(i,j)(1./(alpha+distance(P(i,:),P(j,:))));
fkern = @(i,j)(1./(alpha+pdist2(P(i,:),P(j,:),'euclidean')));

% create the matrix 
K = fkern(1:N,1:N);

% create the H-matrix structure
H = hmtx_cluster(P,'eta',eta,'Nmin',Nmin);

% plot the empty H-matrix
hmtx_plot(H);

% populate H
H = hmtx_fill(H,tol,fkern);

[nB,nE] = hmtx_memory(H);
fprintf('Memory compression: %f\n',nB/getfield(whos('K'),'bytes'));
fprintf('Entry compression: %f\n',nE/(H.nrow*H.ncol));

% plot the H-matrix with the reduced ranks
hmtx_plot(H);
title('H matrix')

% cluster tree
[qr,qc] = hmtx_getclustertree(H);
a = rand(max(qr),3);
figure, scatter3(P(:,1),P(:,2),P(:,3),80,a(qr,:),'filled');

% create fullmatrix from H-matrix (only for small matrices)
Hf = hmtx_full(H);

%% same process with another H-matrix

%Q = [2*rand(N,1)-1 zeros(N,2)];

% kernel function: 1/r, protecting the diagonal with alpha
beta=3;
gkern = @(i,j)(1./(beta+pdist2(P(i,:),P(j,:), 'euclidean').^2));

L = gkern(1:N,1:N);
B = hmtx_cluster(P,'eta',eta,'Nmin',Nmin);
B = hmtx_fill(B,tol,gkern);
hmtx_plot(B);
title('B matrix')


%% check the Frobenius norm of the two matrices
fprintf('Recon. accuracy: ')
relnorm = norm(K-Hf,'fro')/norm(K,'fro');
if relnorm < tol
    fprintf('\ttest passed\n');
else
    fprintf('\t\ttest failed: %e > %e\n',tol);
end

%% transposition
HfT = transpose(Hf);
HT =  test_transpose(H, HfT, tol);

%% matrix inverse
H0 =inv(Hf);
Hinv = test_inv(H,H0,tol);

%% check matrix multiplication and power matrix H*H
H2f =Hf*Hf;
HH = test_mult(H, H2f, tol);

%% addition and subtraction
HpH = test_add(H, Hf, tol);

%% test solver and LU decomposition

%[xL,xU,XR,XL, L, U] = test_LU(H,B);


