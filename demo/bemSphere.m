clear variables, close all

% H matrix params
eta = 2;
Nmin = 20;
tol = 1e-4;
hSize = 0.2;


% surface mesh of a sphere
[P,F] = triSphere(hSize);
M = 2*ones(size(F,1),1);

% data structure
primal = createPrimal2d(P,F,M);
plotScalar2d(primal,'mesh','y');

% electric data
r0 = 1;
V0 = 1;

% the problem has normal derivative as unknowns
% G*du = H*u

% imposed potential on the surface
u = V0*ones(size(primal.Fac2Nod,1),1);

% BEM full matrices
tic
[G0,H0] = bemMatrix2d(primal);
toc

% barycenters
B = barycenter(primal.Node,primal.Fac2Nod);

% G kernel function (from full matrix)
Gkern = @(irow,jcol)G0(irow,jcol);
Hkern = @(irow,jcol)H0(irow,jcol);

% 
tic
G = hmtx_cluster(B,'eta',eta,'Nmin',Nmin);
G = hmtx_fill(G,tol,Gkern);
toc
tic
H = hmtx_cluster(B,'eta',eta,'Nmin',Nmin);
H = hmtx_fill(H,tol,Hkern);
toc

% plot and stats
hmtx_plot(G);
[nB,nE] = hmtx_memory(G);
fprintf('Memory compression: %f\n',nB/(8*G.nrow*G.ncol));
fprintf('Entry compression: %f\n',nE/(G.nrow*G.ncol));

% rhs
b = hmtx_HxM(H,u);

% solution: du = G\b;
errMAX = 1e-10;
iterMAX = 50;

AxFun = @(x)hmtx_HxM(G,x);

% approximate inverse for preconditioning
tol2 = .1; %tol*10^ceil(-log10(tol))*10^(-.5*ceil(-log10(tol)));
Ga = hmtx_coarsening(G,tol2,'kmax',3);
% tic
% invGa = hmtx_inv(Ga);
% toc
% hmtx_plot(invGa);

tic
[L,U] = hmtx_lu(Ga);
toc
LxFun = @(x)hmtx_lsolve(L,x);
UxFun = @(x)hmtx_usolve(U,x);

% gmres w/o preconditioning
tic
[du0,flag0,relres0,iter0,resvec0] = gmres(AxFun,b,[],errMAX,iterMAX,[]);
toc

% gmres w/ preconditioning
tic
[du,flag,relres,iter,resvec] = gmres(AxFun,b,[],errMAX,iterMAX,LxFun,UxFun);
toc

% convergence
h = figure;
hold on
iterMax = max([iter0(2),iter(2)]);
plot(1:iter0(2)+1,log10(relres0*resvec0/resvec0(end)));
plot(1:iter(2)+1,log10(relres*resvec/resvec(end)));
plot(1:iterMax+1,log10(errMAX)*ones(1,iterMax+1),'k--'); % convergence tolerance
axis([0 iterMax+1 log10(errMAX/10) 1])
xlabel('iteration #')
ylabel('log10(|| Ax-b ||/|| b ||)');
legend('w/o precond.','w/ precond');
figureSetup(h);
drawnow

% inverse
tic
invG = hmtx_inv(G);
toc
du = hmtx_HxM(invG,b);

% solve with factorization
tic
[L,U] = hmtx_lu(G);
toc
tic
y = hmtx_lsolve(L,b);
du = hmtx_usolve(U,y);
toc

tic
x = G0\b;
toc

% color maps of solution
plotScalar2d(primal,'field',u,'mesh','y','cLim',[0 1],'title','imposed potential: u');
plotScalar2d(primal,'field',du,'mesh','y','cLim',[0 1],'title','solution: du');

% potential along a line
s = linspace(1.01,10,100)';
Qglo = [s 0*s 0*s];
[G0,H0] = bemMatrix2d(primal,'point',Qglo);
uQ1 = G0*du-H0*u;

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

