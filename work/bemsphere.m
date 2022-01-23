clear variables, close all

% surface mesh of a sphere
Ps = triSphere(.1);

% tetrahedralization
DT = DelaunayTri([Ps; 0 0 0]);
P = DT.X;
T = DT.Triangulation;
M = 2*ones(size(T,1),1);

% data structure
primal = createPrimal3d(P,T,M);
plotScalar3d(primal,'mesh','y');

% electric data
r0 = 1;
V0 = 1;

% the problem has normal derivative as unknowns
% G*du = H*u

% imposed potential on the surface
fsurf = findsurf3d(primal);
usurf = ones(length(fsurf),1);

% BEM matrices
tic
[G,H] = bem_matrix3d(primal);
toc

% rhs
b = H*usurf;

% % calculate G with function fkern
B = barycenter(primal.node,primal.face_node(fsurf,:));
fkern = @(irow,jcol)bem_Gcoeff(irow,jcol,primal.node,primal.face_node(fsurf,:),B);
% 
tic
Gk = hmtx_cluster(B,'eta',2,'Nmin',50);
toc
% hmtx_plot(Gk);
eps = 1e-4;
tic
Gk = hmtx_fill(Gk,eps,fkern);
toc
% tic
% Hk = hmtx_cluster_class(B,'eta',2,'Nmin',50);
% Hk = hmtx_fill(Hk,eps,fkern);
% toc

return

% solution
% dusurf = G\b;
errMAX = 1e-6;
M = diag(diag(G));
%[dusurf,flag,relres,iter,resvec] = gmres(G,b,[],errMAX,[]);
afun = @(x)Axfun(x,Gk);
[dusurf,flag,relres,iter,resvec] = gmres(afun,b,[],errMAX,[]);

% convergence
h = figure;
plot(1:iter(2)+1,log10(relres*resvec/resvec(end)),'r-');
hold on
plot(1:iter(2)+1,log10(errMAX)*ones(1,iter(2)+1),'b--'); % convergence tolerance
axis([0 iter(2)+1 log10(errMAX/10) 1])
xlabel('iteration #')
ylabel('log10(|| Ax-b ||/|| b ||)');
figuresetup(h);

% set white background
set(gcf,'Color',[1 1 1]);

% preallocate
u = zeros(primal.face_num,1);
du = zeros(primal.face_num,1);

% solution projected to volumic data structure
u(fsurf) = usurf;
du(fsurf) = dusurf;

% color maps of solution
plotscalar3d(primal,'field',u,'mesh','y','cLim',[0 1],'title','solution: u');
plotscalar3d(primal,'field',du,'mesh','y','cLim',[0 1],'solution: du');

% potential along a line
s = linspace(1.01,10,100)';
Qglo = [s 0*s 0*s];
[G,H] = bem_matrix3d(primal,'point',Qglo);
fsurf = findsurf3d(primal);
uQ1 = G*du(fsurf)-H*u(fsurf);

% electric field E=-grad(u)
E = bem_vectfield3d(primal,du,Qglo);

r = sqrt(sum(Qglo.^2,2));
h1 = figure;
subplot(2,1,1), plot(r,V0*r0./r,'r');
hold on, plot(r,uQ1,'g');
xlabel('distance (m)'); ylabel('potential (V)'); legend('analytic','bem');
figuresetup(h1)

subplot(2,1,2), plot(r,V0*r0./r.^2,'r');
hold on, plot(r,E(:,1),'g');
hold on, plot(r,-interpderiv1(r,uQ1,r),'b');
xlabel('distance (m)'); ylabel('electric field (V/m)'); legend('analytic','bem','numeric gradient');
figuresetup(h1)

