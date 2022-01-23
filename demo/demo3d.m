clear variables, close all

% matrix generation
alpha = 1e-2;
N = 500;
tol = 1e-4;  % tolerance of the ACA method
eta = 3;
Nmin = 10;

% create a grid of N elements in [-1,1]. The clustering algorithm works
% with 3d points
P = bsxfun(@times,rand(N,3),[1 3 5]);

figure,plot3(P(:,1),P(:,2),P(:,3),'o');

% kernel function: 1/r, protecting the diagonal with alpha
fkern = @(i,j)(1./(alpha+distance(P(i,:),P(j,:))));

% create the matrix 
K = fkern(1:N,1:N);

% create the H-matrix structure
H = hmtx_cluster(P,'eta',eta,'Nmin',Nmin);

% plot the empty H-matrix
hmtx_plot(H);

% populate H
H = hmtx_fill(H,tol,fkern);

% plot the H-matrix with the reduced ranks
hmtx_plot(H);
%print('-depsc2','-r600','H.eps');
% create fullmatrix from H-matrix (only for small matrices)
Hf = hmtx_full(H);
H0 = fkern(1:N,1:N);

%% check the Frobenius norm of the two matrices
fprintf('Recon. accuracy: ')
relnorm = norm(K-Hf.M,'fro')/norm(K,'fro');
if relnorm < tol
    fprintf('\ttest passed\n');
else
    fprintf('\t\ttest failed: %f > %f\n',tol);
end

%% check matrix multiplication
% generate a random matrix
R = rand(N);

% left product R*H
fprintf('Left product: ')
relnorm = norm(R*Hf.M-hmtx_mtimes(R,H),'fro')/norm(R*Hf.M,'fro');
if relnorm < tol
    fprintf('\t\ttest passed\n');
else
    fprintf('\t\ttest failed: %f > %f\n',relnorm,tol);
end

% right product H*R
fprintf('Right product: ')
relnorm = norm(Hf.M*R-hmtx_mtimes(H,R),'fro')/norm(Hf.M*R,'fro');
if relnorm < tol
    fprintf('\t\ttest passed\n');
else
    fprintf('\t\ttest failed: %f > %f\n',relnorm,tol);
end

% % power matrix H*H
% fprintf('Self product: ')
% H2 = hmtx_full(hmtx_mtimes(H,H));
% relnorm = norm(Hf.M*Hf.M-H2.M,'fro')/norm(Hf.M*Hf.M,'fro');
% if relnorm < 100*eps
%     fprintf('\t\ttest passed\n');
% else
%     fprintf('\t\ttest failed: %f > %f\n',relnorm,100*eps);
% end


%% addition and subtraction
fprintf('Addition: ')
HpH = hmtx_full(hmtx_add(H,H,'+'));
relnorm = norm(Hf.M+Hf.M-HpH.M,'fro')/norm(Hf.M+Hf.M,'fro');
if relnorm < tol
    fprintf('\t\ttest passed\n');
else
    fprintf('\t\ttest failed: %f > %f\n',relnorm,tol);
end

fprintf('Subtraction: ')
HmH = hmtx_full(hmtx_add(H,H,'-'));
relnorm = norm(HmH.M,'fro');
if relnorm < tol
    fprintf('\t\ttest passed\n');
else
    fprintf('\t\ttest failed: %f > %f\n',relnorm,tol);
end

%% LU

test_LU(H, [], tol)

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


% %% test new Lsolve
% beta=3;
% gkern = @(i,j)(1./(beta+pdist2(P(i,:),P(j,:), 'euclidean')*5));
% 
% Li = gkern(1:N,1:N);
% B = hmtx_cluster(P,'eta',eta,'Nmin',Nmin);
% B = hmtx_fill(B,tol,gkern);
% hmtx_plot(B);
% 
% title('B matrix')
% 
% test_newsolve(H, H0,B, tol);


%% coarsening matrix tests
newtol = tol * 10^3;
test_coarse(H, H0, tol, newtol);

return


