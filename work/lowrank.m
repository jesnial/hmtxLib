clear variables, close all

% number of points per direction
N = 1000;
% parameter for non-singularity
alpha = 1e-1;

x = linspace(-1,1,N);
y = linspace(-1,1,N);
[X,Y] = meshgrid(x,y);

% kernel matrix
K = 1./(alpha+abs(X-Y));

% singular value decomposition
[U,S,V] = svd(K);
sigma = diag(S);

% low rank reconstruction
err = zeros(N,1);
for r = 1:N
    Kr = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    err(r) = norm(K-Kr,'fro')/norm(K,'fro');
end

% plot
h = figure;
subplot(1,2,1), semilogy(sigma);
subplot(1,2,2), semilogy(err);

