clear variables, close all

%% parameters
% cylinder: nuber of circumference divisions
n1 = 30;
% cylinder: number of height divisions
nz = 10;
% cylinder: radius
r1 = .175;
% cylinder: semi-height
zmax = .193;
% coil: radius
r2 = .1255;
% coil: nuber of circumference divisions
n2 = 20;
% frequency
f = 1000;
% current
I0 = 2;

eps = 0.5;  % tolerance of the ACA method
eta = 1e-4;    % threshold for admissibility
Nmin = 32;   % minimum number of points in a leaf of the cluster tree
tol = 1e-8;  % tolerance for GMRES

%% create cylinder
[P,T] = thincylinder(r1,-zmax,zmax,n1,nz);
M = ones(size(T,1),1);
 
%% create coil
% angle
ad = linspace(0,2*pi,n2);
% points
Q = [r2*cos(ad(:)) r2*sin(ad(:)) 0*ad(:)];
% start points
coil.P1 = Q(1:end-1,:);
% end points
coil.P2 = Q(2:end,:);
coil.I0 = I0;

%% tracks
track = [];

%% create data structures
primal = createprimal2d(P,T,M);
primal.scale = 1;
primal.face_node = int32(primal.face_node);
primal.sort = int32(ones(primal.face_num,1));
mat = setmatproperties(M);
clear P T M

h = plotscalar2d(primal,'mesh','y');
plotline3d(coil.P1,coil.P2,'handle',h,'arrow','y'); drawnow

%% material properties
mat.sigma(1) = 59.5e6;
mat.delta(1) = .5e-3;

%% define boundary conditions
[idfix,xfix,idfloat] = eddyshellbc(primal,track);

%% calculate coefficients
% calculate inductance matrix
L = inductancematrix(primal,mat);
% calculate resistance
R = resistancematrix(primal,mat);
% coil/structure mutal coupling
b = coilshield(primal,coil);

%% harmonic solution
x = eddyshellsolve(R,L,b,f,idfix,xfix,idfloat);

%% ACA-Hmatrix approach
clear mim_Lcoeff % necessary to erease the persistent variable inside the funciton
fkern1 = @(irow,jcol)mim_Lcoeff(irow,jcol,primal);
fkern2 = @(irow,jcol)mim_Lcoeff_num(irow,jcol,primal.node,primal.face_node);


% filling matrix
tic; Lk = hmtx_cluster(primal.node,'eta',eta,'Nmin',Nmin); toc
tic; Lk = hmtx_fill(Lk,eps,fkern1,fkern2); toc
x2 = eddyshellsolve_Hmatrix(R,Lk,b,f,idfix,xfix,idfloat,'gout','n');
% hmtx_plot(Lk);

