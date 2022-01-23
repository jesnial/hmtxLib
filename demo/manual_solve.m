clear variables, close all

%% parameters
alpha = 1e-2;  % singularity shift
N = 8;        % number of points
R0 = 1;
tol = 1;    % tolerance of the ACA method
eta = 3;       % admissibility criterion
Nmin = 4;      % minimum points per cluster

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

%% create triangular matrix

Hf = getfield(hmtx_full(H), 'M');

% create L and U matrices
L = hmtx_tril(H);
U = hmtx_triu(H);

hmtx_plot(L);
title('L matrix')
hmtx_plot(U);
title('U matrix')

Lf = getfield(hmtx_full(L),'M');
Uf = getfield(hmtx_full(U),'M');

%% manually solve Lx = b

b = ones(L.nrow,1);
x1 = hmtx_lsolve(L,b);
xf = Lf\b;

% global supermatrix
x = zeros(L.ncol,1);

% upper block
x(L.M{1,1}.M{1,1}.jcol,:) = L.M{1,1}.M{1,1}.M\b(L.M{1,1}.M{1,1}.irow,:);
L1121x = zeros(L.M{1,1}.M{2,1}.nrow,1);
L1121x(L.M{1,1}.M{2,1}.irow,:) =  hmtx_HxM(L.M{1,1}.M{2,1},x);
%L.M{1,1}.M{2,1}.U*L.M{1,1}.M{2,1}.V'
%L1121x(L.M{1,1}.M{2,1}.irow,:) =  getfield(hmtx_full(L.M{1,1}.M{2,1}),'M')*x(L.M{1,1}.M{2,1}.jcol,:);
x(L.M{1,1}.M{2,2}.jcol,:) = L.M{1,1}.M{2,2}.M\(b(L.M{1,1}.M{1,1}.irow,:)-L1121x(L.M{1,1}.M{1,1}.irow,:));
L1121x
x-xf

% middle block solution
% rhs
L21x = zeros(L.M{2,1}.nrow,1);
L21x(L.M{2,1}.irow,:) =  hmtx_HxM(L.M{2,1},x);
%L.M{2,1}.U*L.M{2,1}.V'

%L21x(L.M{2,1}.irow,:) =  getfield(hmtx_full(L.M{2,1}),'M')*x(L.M{2,1}.jcol,:);
L21x
x-xf

% lower block solution
b2 = b -L21x;
x(L.M{2,2}.M{1,1}.jcol,:) = L.M{2,2}.M{1,1}.M\b2(L.M{2,2}.M{1,1}.irow,:);
L2221x = zeros(L.M{2,2}.M{2,1}.nrow,1);
L2221x(L.M{2,2}.M{2,1}.irow,:) =  hmtx_HxM(L.M{2,2}.M{2,1},x);
%L.M{2,2}.M{2,1}.U*L.M{2,2}.M{2,1}.V'
%L2221x(L.M{2,2}.M{2,1}.irow,:) =  getfield(hmtx_full(L.M{2,2}.M{2,1}),'M')*x(L.M{2,2}.M{2,1}.jcol,:);
x(L.M{2,2}.M{2,2}.jcol,:) = L.M{2,2}.M{2,2}.M\(b2(L.M{2,2}.M{1,1}.irow,:)-L2221x(L.M{2,2}.M{1,1}.irow,:));
L2221x
x-xf

%% comparisons
norm(x-xf)/norm(xf)
norm(x-x1)/norm(x1)
norm(x1-xf)/norm(xf)

