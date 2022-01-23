%% Spinterometer1: solution with explicit calculation of inv(Gk)
% BEM formulation of two spheres with imposed potentials 'u'
% the problem has normal derivative of potential 'du' as unknowns
% G*du = H*u = b

clear variables, close all

%% parameters
r0 = 1.5;    % sphere radius
V0 = 2;      % voltage applied between the spheres
d = 4;       % distance between the sphere centers
eps = 1e-6;  % tolerance of the ACA method

%% data structure: calculation of nodes P, faces F, barycenters B
% surface mesh of a sphere
[P,F] = trisphere(.4);

% normalize radius to r0
r = sqrt(sum(P.^2,2));
P = r0*P./repmat(r,1,3);

% add second sphere
F = [F; F+size(P,1)];
P = [P; [P(:,1)+d P(:,2:3)]];

% barycenters of triangles
B = barycenter(P,F);

%% matrix calculation
% calculate G with function fkern
fkern1 = @(irow,jcol)bem_Gcoeff(irow,jcol,P,F,B);
fkern2 = @(irow,jcol)bem_Hcoeff(irow,jcol,P,F,B);

tic; Gk = hmtx_cluster(B,'eta',1.5,'Nmin',50); toc
Hk = Gk;

% populate Gk
tic; Gk = hmtx_fill(Gk,eps,fkern1); toc

% populate Hk
tic; Hk = hmtx_fill(Hk,eps,fkern2);
% fix diagonal entries
Hk = hmtx_fixdiag(Hk); toc

% rhs: imposed potential on the surface
u = V0/2*ones(size(F,1),1);
u(B(:,1) > d/2) = -V0/2;

%% solution du = inv(Gk)*b
% matrix inverse
tic; Gkinv = hmtx_inv(Gk); toc

% rhs
b = hmtx_HxM(Hk,u);

% solution
du = hmtx_HxM(Gkinv,b);

% color maps of solution
figure; patch('Faces',F,'Vertices',P,'FaceVertexCData',du,'FaceColor','flat','EdgeColor','k');
view([1 1 1]); axis equal; colorbar