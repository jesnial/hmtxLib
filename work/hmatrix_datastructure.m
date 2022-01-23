clear variables, close all
M = 4000;
N = 4000;

% vector
x = rand(N,1);

% matrix
A = cell(2,2);

A{1,1}.type = 'fullmatrix';
A{1,1}.m = M/2;
A{1,1}.n = N/2;
A{1,1}.M = rand(A{1,1}.m,A{1,1}.n);

A{1,2}.type = 'fullmatrix';
A{1,2}.m = M/2;
A{1,2}.n = N/2;
A{1,2}.M = rand(A{1,2}.m,A{1,2}.n);

A{2,1}.type = 'fullmatrix';
A{2,1}.m = M/2;
A{2,1}.n = N/2;
A{2,1}.M = rand(A{2,1}.m,A{2,1}.n);

% block ACA
% radius of sources
R0 = 1;
% inner radius of targets
Rin = 2;
% outer radius of targets
Rout = 4;
% source points
r = R0*rand(N/2,1);
theta = pi*rand(N/2,1);
phi = 2*pi*rand(N/2,1);
P = [r.*sin(theta).*cos(phi) r.*sin(theta).*sin(phi) r.*cos(theta)];
% target points
r = Rin+(Rout-Rin)*rand(M/2,1);
theta = pi*rand(M/2,1);
phi = 2*pi*rand(M/2,1);
Q = [r.*sin(theta).*cos(phi) r.*sin(theta).*sin(phi) r.*cos(theta)];

A{2,2}.type = 'rkmatrix';
A{2,2}.m = M/2;
A{2,2}.n = N/2;
K = 1./distance(Q,P);
A{2,2}.M = K;

R = K;
normK = norm(K,'fro')^2;
err = 1;
U = [];
V = [];
eps = 1e-6;
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
    err = normR/normK;
end
A{2,2}.r = size(U,1);
A{2,2}.U = U;
A{2,2}.VT = V;

% assemble full matrix
Am = [A{1,1}.M A{1,2}.M; A{2,1}.M A{2,2}.M];
tic
b = Am*x;
toc

tic
b2 = hmtx_Ax(A,x);
toc

