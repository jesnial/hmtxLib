clear variables, close all

% Alg 3.1, S. Rjasanow, O. Steinbach, "The fast solution of Boundary Integral
% Equations", Springer, 2007

% number of points per direction
N = 5;

y = linspace(-1,-.5,N);
x = linspace(.5,1,N);
[X,Y] = meshgrid(x,y);

% kernel matrix
K = 1./(abs(X-Y));

% initialization
R = K;
S = zeros(size(K));
err = 1;

eps = 1e-2;
normK = norm(K,'fro');

U = [];
V = [];
while err > eps
    % find pivot
    [maxK,imaxR] = max(R(:));
    [imax,jmax] = ind2sub(size(K),imaxR);
    % normalizing constant
    gamma = 1./R(imaxR);
    % new functions
    u = gamma*R(:,jmax);
    v = R(imax,:);
    % new residual
    R = R-u*v;
    % new approximation
    U = [U u];
    V = [V; v];
    normR = norm(R,'fro');
    err = normR/normK
end

[U2,VT2,err2] = aca_fullpivot(K,eps)

