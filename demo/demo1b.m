clear variables, close all

%% parameters
alpha = 1e-2;  % singularity shift
N = 1200;        % number of points
R0 = 1;
tol = 1e-6;    % tolerance of the ACA method
eta = 3;       % admissibility criterion
Nmin = 20;      % minimum points per cluster

% theta = linspace(-pi/2,pi/2,round(N/3));
% phi = linspace(0,2*pi,round(2*N/3));
% [R,THETA,PHI] = meshgrid(R0,theta,phi);
% P = zeros(numel(R),3);
% [P(:,1),P(:,2),P(:,3)] = sph2cart(R(:),THETA(:),PHI(:));
% P = triSphere(0.2);

% create a grid of N elements in [-1,1] (in 3d for the clustering algorithm)
P = [logspace(-1,1,N).' zeros(N,2)];
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


%% test new Lsolve
beta=3;
gkern = @(i,j)(1./(beta+pdist2(P(i,:),P(j,:), 'euclidean')*5));

Li = gkern(1:N,1:N);
B = hmtx_cluster(P,'eta',eta,'Nmin',Nmin);
B = hmtx_fill(B,tol,gkern);
hmtx_plot(B);

title('B matrix')

test_newsolve(H, H0,B, tol);


%% test solver and LU decomposition

test_LU(H, H0, tol);

%% coarsening matrix tests
newtol = tol * 10^3;
test_coarse(H, H0, tol, newtol);

return
%% test cholesky
HH = hmtx_herm(H);
HHH = hmtx_mult(H,HH , H);

hmtx_issym(HHH, tol)

HHHf = getfield(hmtx_full(HHH), 'M');
HHf = Hf';
HH_f = getfield(hmtx_full(HH), 'M');
fprintf('herm: ')
relnorm = norm(HHf-HH_f,'fro')/norm(HHf,'fro');
if relnorm < tol
    fprintf('\t\t\t\ttest passed: %e < %e',relnorm,tol);
else
    fprintf('\t\t\t\ttest failed: %e > %e',relnorm,tol);
end

tic;
Lh = hmtx_cholesky(HHH);
t0 = toc;

[L, flag] = chol(HHHf, 'lower');
flag

Lh_f = getfield(hmtx_full(Lh), 'M');

figure(20);
imshow(L);
title('full L')

hmtx_plot(Lh);
title('H-matrix L')

figure(21);
imshow(Lh_f);
title('hmtx to full L')


fprintf('Cholesky: ')
relnorm = norm(L*L'-Lh_f*Lh_f','fro')/norm(L*L','fro');
if relnorm < tol
    fprintf('\t\t\t\ttest passed: %e < %e',relnorm,tol);
else
    fprintf('\t\t\t\ttest failed: %e > %e',relnorm,tol);
end
fprintf('\ttime: %f s\n',t0);



    