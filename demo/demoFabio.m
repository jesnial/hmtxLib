clear all, close all

%% parameters
alpha = .1;    % singularity shift
tol = 1e-4;    % tolerance of the ACA method
kMax = Inf;
eta = 4;       % admissibility criterion
N = 80;
Np = N*(N-1)+2; % number of actual points
Nmin = min([32,max([round(Np/10),3])]);      % minimum points per cluster

[X,Y,Z] = sphere(N);
P = uniquetol([X(:) Y(:) Z(:)],'ByRows',true);

% rng(1);
% P = rand(N^2,3);

N = size(P,1);        % number of points


% create a grid of N elements in [-1,1] (in 3d for the clustering algorithm)
% P = [linspace(-1,1,N).' zeros(N,2)];
% P = [2*rand(N,1)-1 zeros(N,2)];
% P = rand(N,3);

% kernel function: 1/r, protecting the diagonal with alpha
fKern = @(i,j)1./(alpha+pdist2(P(i,:),P(j,:),'euclidean')).^2;

% create the full matrix 
H0 = fKern(1:N,1:N);

% create the H-matrix structure
H = hmtx_cluster(P,'eta',eta,'Nmin',Nmin);

% plot the H-matrix with the reduced ranks
hmtx_plot(H);
title('H-matrix structure')

% populate H
H = hmtx_fill(H,tol,fKern,'kMax',kMax);
H = hmtx_tH(.5,hmtx_add(H,hmtx_transpose(H),'+'));

% plot the H-matrix with the reduced ranks
hmtx_plot(H);
title('H-matrix')

% memory data
[nB,nE] = hmtx_memory(H);
fprintf('Matrix size: %d x %d\n',H.nrow,H.ncol);
fprintf('Memory compression: %f\n',nB/(8*H.nrow*H.ncol));
fprintf('Entry compression: %f\n',nE/(H.nrow*H.ncol));

% cluster tree
[qr,qc] = hmtx_getclustertree(H);
a = rand(max(qr),3);
figure
scatter3(P(:,1),P(:,2),P(:,3),80,a(qr,:),'filled');
axis equal;

% create fullmatrix from H-matrix (only for small matrices)
Hf = getfield(hmtx_full(H),'M');

%% check the Frobenius norm of the two matrices
relnorm = norm(H0-Hf,'fro')/norm(H0,'fro');
fprintf('Reconstruction accuracy: ')
if relnorm < tol
    fprintf('\ttest passed: %e < %e\n',relnorm,tol);
else
    fprintf('\ttest failed: %e > %e\n',relnorm,tol);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % use H as rhs
% b = getfield(hmtx_full(H),'M');
% 
% % extract trilower matrix
% L = hmtx_tril(H);
% 
% % convert to full
% Lf = getfield(hmtx_full(L),'M');
% 
% % solution with full matrix
% xf = Lf\b;
% 
% % H-matrix solve
% 
% tic;
% X = hmtx_leftsolve(L,H);
% t0 = toc;
% 
% % check 
% fprintf('left-solve: ')
% relnorm = norm(getfield(hmtx_full(X),'M')-xf,'fro')/norm(xf,'fro');
% if relnorm < tol
%     fprintf('\t\t\ttest passed: %e < %e',relnorm,tol);
% else
%     fprintf('\t\t\ttest failed: %e > %e',relnorm,tol);
% end
% fprintf('\ttime: %f s\n',t0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % extract triupper matrix
% U = hmtx_triu(H);
% 
% % convert to full
% Uf = getfield(hmtx_full(U),'M');
% 
% % solution with full matrix
% xUf = b/Uf;
% 
% % H-matrix solve
% 
% tic;
% XU = hmtx_rightsolve(U,H);
% t0 = toc;

% check 
% fprintf('right-solve: ')
% relnorm = norm(getfield(hmtx_full(XU),'M')-xUf,'fro')/norm(xUf,'fro');
% if relnorm < tol
%     fprintf('\t\t\ttest passed: %e < %e',relnorm,tol);
% else
%     fprintf('\t\t\ttest failed: %e > %e',relnorm,tol);
% end
% fprintf('\ttime: %f s\n',t0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% transposition
% test_transpose(H,H0,tol);
% 
%% scalar multiplication
% test_scalar(H,H0,tol);

%% addition and subtraction
 test_add(H,H0,tol);
% test_sub(H,H0,tol);

%% combination: make matrix symmetric
% tic;
% Hs = hmtx_tH(.5,hmtx_add(H,hmtx_transpose(H),'+'));
% t0 = toc;
% 
% H0s = (H0+H0')/2;
% relnorm = norm(getfield(hmtx_full(Hs),'M')-H0s,'fro')/norm(H0,'fro');
% fprintf('Symmetrization: ')
% if relnorm < tol
%     fprintf('\t\ttest passed: %e < %e',relnorm,tol);
% else
%     fprintf('\t\ttest failed: %e > %e',relnorm,tol);
% end
% fprintf('\ttime: %f s\n',t0);
% 
% %% check matrix multiplication
% test_matvect(H,H0,tol);

% %% l-solve
% test_lusolve(H,H0,tol);
% 
% %% check H-matrix multiplication
% test_mult(H,H0,tol);

% %% check inverse H-matrix
% test_inv(H,H0,tol);


% %% test LU
% test_LU(H, H0, tol);

%% test chol
% tic
% [L0,U0] = hmtx_lu(H);
% toc
% % tic
% % L = hmtx_chol(H);
% % toc
% b = ones(H.nrow,1);
% tic; x0 = hmtx_usolve(U0,hmtx_lsolve(L0,b)); toc;
% tic; x = hmtx_ltsolve(L,hmtx_lsolve(L,b)); toc;
% tic; x1 = hmtx_usolve(hmtx_ctranspose(L),hmtx_lsolve(L,b)); toc;
% 
% figure,plot(x,'o'),hold on,plot(x0,'*')
% return
% %% coarsening 
% newtol = tol * 10^3;
% test_coarse(H, H0, tol, newtol);

