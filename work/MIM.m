clear variables, close all

%% parameters
r0 = 1.5;    % sphere radius
V0 = 2;      % voltage applied between the spheres
d = 4;       % distance between the sphere centers
eps = 1e-2;  % tolerance of the ACA method
tol = 1e-10; % tolerance for GMRES

%% data structure: calculation of nodes P, faces F, barycenters B
% surface mesh of a sphere
[P,F] = trisphere(.2);

% normalize radius to r0
r = sqrt(sum(P.^2,2));
P = r0*P./repmat(r,1,3);

% add second sphere
F = [F; F+size(P,1)];
P = [P; [P(:,1)+d P(:,2:3)]];

r0 = 1.5;    % track radius
nt = 60;     % nuber of circumference divisions
ntz = 25;     % number of height divisions
zmax = .2;  % track: semi-height

%% mesh data structure
% create tracks
[P,F] = thincylinder(r0,-zmax,zmax,nt,ntz);
M = [ones(size(F,1),1)];



% barycenters of triangles
B = barycenter(P,F);

%% matrix calculation
% calculate G with function fkern
fkern1 = @(irow,jcol)mim_Lcoeff(irow,jcol,P,F);
fkern2 = @(irow,jcol)mim_Lcoeff_num(irow,jcol,P,F);

tic; Gk = hmtx_cluster(P,'eta',2,'Nmin',32); toc
hmtx_plot(Gk);

% populate Gk
tic; Gk = hmtx_fill(Gk,eps,fkern1); toc
hmtx_plot(Gk)

% create corresponding full matrix
L2 = hmtx_full(Gk);

irow = 1:size(P,1);
jcol = 1:size(P,1);
tic; L3 = newinductancematrixf(P,int32(F),int32(irow(:)),int32(jcol(:))); toc

tic; Lf = inductancematrixf(P,int32(F)); toc

norm(L2.M-Lf)/norm(Lf)
