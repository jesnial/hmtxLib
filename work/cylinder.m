%% cylinder: GMRES solution
% BEM formulation of two cylindrical electrodes with imposed potentials 'u'
% the problem has normal derivative of potential 'du' as unknowns
% G*du = H*u = b

clear variables, close all

%% parameters
r1 = .5;     % inner radius
r2 = .7;     % outer radius
dz = 2;      % heigth of cylinder
V0 = 2;      % voltage applied between the eletrodes
eps = 1e-4;  % tolerance of the ACA method
tol = 1e-10; % tolerance for GMRES
h = 0.07;    % mesh size;

%% data structure: calculation of nodes P, faces F, barycenters B
% number of divisions
na1 = round(2*pi*r1/h);
na2 = round(2*pi*r2/h);
nz = round(dz/h);

% surface mesh of inner cylinder
[P1,F1] = thincylinder(r1,-dz/2,dz/2,na1,nz);

% surface mesh of outer cylinder
[P2,F2] = thincylinder(r2,-dz/2,dz/2,na2,nz);

% concatenate structure
F = [F1; F2+size(P1,1)];
P = [P1; P2];

% barycenters of triangles
B = barycenter(P,F);

%% matrix calculation
% calculate G with function fkern
fkern1 = @(irow,jcol)bem_Gcoeff(irow,jcol,P,F,B);
fkern2 = @(irow,jcol)bem_Hcoeff(irow,jcol,P,F,B);

tic; Gk = hmtx_cluster(B,'eta',3,'Nmin',50); toc
Hk = Gk;

% populate Gk
tic; Gk = hmtx_fill(Gk,eps,fkern1); toc

% populate Hk
tic; Hk = hmtx_fill(Hk,eps,fkern2);
% fix diagonal entries
Hk = hmtx_fixdiag(Hk); toc

% rhs: imposed potential on the surface
u = V0/2*ones(size(F,1),1);
idx = sqrt(B(:,1).^2+B(:,2).^2) < (r1+r2)/2;
u(idx) = -V0/2;

%% solution with GMRES
% rhs
tic; b = hmtx_HxM(Hk,u); toc;

% anonymous function to calculate A*x
afun = @(x)hmtx_HxM(Gk,x);
tic; [du,flag,relres,iter,resvec] = gmres(afun,b,[],1e-10,100); toc;

% plot GMRES convergence
figure; plot(0:iter(2),log10(resvec/norm(b)),'-o');
xlabel('Iteration number');
ylabel('Relative residual');

% color maps of solution
figure; patch('Faces',F,'Vertices',P,'FaceVertexCData',du,'FaceColor','flat','EdgeColor','k');
view([1 1 1]); axis equal; colorbar