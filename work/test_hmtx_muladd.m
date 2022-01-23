clear variables, close all

%% parameters
r0 = 1.5;    % sphere radius
V0 = 2;      % voltage applied between the spheres
eps = 1e-4;  % tolerance of the ACA method
tol = 1e-10; % tolerance for GMRES
dx = 2.1*r0;
dy = 3*r0;

%% data structure: calculation of nodes P, faces F, barycenters B
% surface mesh of a sphere
[P,F] = trisphere(.5);

% normalize radius to r0
r = sqrt(sum(P.^2,2));
P = r0*P./repmat(r,1,3);

% add second sphere
F = [F; F+size(P,1); F+2*size(P,1); F+3*size(P,1);];
P = [
    P(:,1)-dx/2 P(:,2)-dy/2 P(:,3)
    P(:,1)+dx/2 P(:,2)-dy/2 P(:,3)
    P(:,1)-dx/2 P(:,2)+dy/2 P(:,3)
    P(:,1)+dx/2 P(:,2)+dy/2 P(:,3)
    ];

% barycenters of triangles
B = barycenter(P,F);

%% matrix calculation
% calculate G with function fkern
fkern1 = @(irow,jcol)bem_Gcoeff(irow,jcol,P,F,B);

tic; Gk = hmtx_cluster(B,'eta',3,'Nmin',100); toc
tic; Gk = hmtx_fill(Gk,eps,fkern1); toc

% color maps of solution
h2 = figure; patch('Faces',F,'Vertices',P,'FaceColor','r','EdgeColor','k');
view([1 1 1]); axis equal; colorbar
figuresetup(h2);

M = Gk;

% A11 = Gk.M{1,1};
% A12 = Gk.M{1,2};
% A21 = Gk.M{2,1};
% A22 = Gk.M{2,2};
% 
% A11f = hmtx_full(A11);
% A12f = hmtx_full(A12);
% A21f = hmtx_full(A21);
% A22f = hmtx_full(A22);
% 
% M = hmtx_muladd2([],A12,A21);
% Mf = A12f.M*A21f.M;
% norm(Mf-M.M,'fro')

Gf = hmtx_full(M);
tic; D = Gf.M+Gf.M*Gf.M; toc;
tic; Dk = hmtx_muladd2(M,M,M); toc;
F = hmtx_add(M,hmtx_muladd2([],M,M),1,1);
tic; S = Gf.M+Gf.M; toc;
tic; Sk = hmtx_add(M,M,1,1); toc;

hmtx_plot(M);
hmtx_plot(Dk);
Df = hmtx_full(Dk);
Sf = hmtx_full(Sk);
norm(Df.M-D,'fro')
norm(Sf.M-S,'fro')
