clear variables, close all

% H matrix params
eta = 2;
Nmin = 20;
tol = 1e-4;

% surface mesh of a sphere
[P0,F0] = triSphere(.3);
M0 = ones(size(F0,1),1);

% two spheres
P = [bsxfun(@plus,P0,[1.5 0 0]); bsxfun(@plus,P0,[-1.5 0 0])];
F = [F0; F0+size(P0,1)];
M = [2*M0; 3*M0];

% data structure
primal = createPrimal2d(P,F,M);
plotScalar2d(primal,'mesh','y');

% electric data
r0 = 1;
V0 = 1;

% the problem has normal derivative as unknowns
% G*du = H*u

% imposed potential on the surface
u = [V0*ones(size(M0)); -V0*ones(size(M0))];

% BEM full matrices
tic
[G,H] = bemMatrix2d(primal);
toc

% barycenters
B = barycenter(primal.Node,primal.Fac2Nod);

% G kernel function (from full matrix)
Gkern = @(irow,jcol)G(irow,jcol);
Hkern = @(irow,jcol)H(irow,jcol);

% 
tic
Gk = hmtx_cluster(B,'eta',eta,'Nmin',Nmin);
Gk = hmtx_fill(Gk,tol,Gkern);
toc
tic
Hk = hmtx_cluster(B,'eta',eta,'Nmin',Nmin);
Hk = hmtx_fill(Hk,tol,Hkern);
toc

% plot and stats
hmtx_plot(Gk);
[nB,nE] = hmtx_memory(Gk);
fprintf('Memory compression: %f\n',nB/(8*Gk.nrow*Gk.ncol));
fprintf('Entry compression: %f\n',nE/(Gk.nrow*Gk.ncol));

% rhs
b = hmtx_HxM(Hk,u);

% solution: du = G\b;
errMAX = 1e-8;
iterMAX = 50;

AxFun = @(x)hmtx_HxM(Gk,x);

% P = hmtx_makesparse(Gk);
% [L,U] = lu(P);
% PxFun = @(x)hmtx_HxM(invGk,x);

tic
[du,flag,relres,iter,resvec] = gmres(AxFun,b,[],errMAX,iterMAX,PxFun);
toc

% convergence
h = figure;
plot(1:iter(2)+1,log10(relres*resvec/resvec(end)),'r-');
hold on
plot(1:iter(2)+1,log10(errMAX)*ones(1,iter(2)+1),'b--'); % convergence tolerance
axis([0 iter(2)+1 log10(errMAX/10) 1])
xlabel('iteration #')
ylabel('log10(|| Ax-b ||/|| b ||)');
figureSetup(h);

% set white background
set(gcf,'Color',[1 1 1]);

% % inverse
% tic
% invGk = hmtx_inv(Gk);
% toc
% du = hmtx_HxM(invGk,b);
% 
% tic
% x = G\b;
% toc

% color maps of solution
plotScalar2d(primal,'field',u,'mesh','y','title','imposed potential: u');
plotScalar2d(primal,'field',du,'mesh','y','title','solution: du');

return
% potential along a line
s = linspace(1.01,10,100)';
Qglo = [s 0*s 0*s];
[G,H] = bemMatrix2d(primal,'point',Qglo);
uQ1 = G*du-H*u;

% electric field E = -grad(u)
E = -bemVectorField2d(primal,u,du,Qglo);

r = sqrt(sum(Qglo.^2,2));
h1 = figure;
subplot(2,1,1), plot(r,V0*r0./r,'r*');
hold on, plot(r,uQ1,'go');
xlabel('distance (m)'); ylabel('potential (V)'); legend('analytic','bem');
figureSetup(h1)

pDeriv = interpolateDerivative1(r,uQ1);
subplot(2,1,2), plot(r,V0*r0./r.^2,'r*');
hold on, plot(r,E(:,1),'go');
hold on, plot(r,-ppval(pDeriv,r),'bs');
xlabel('distance (m)'); ylabel('electric field (V/m)'); legend('analytic','bem','numeric gradient');
figureSetup(h1)

