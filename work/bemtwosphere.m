clear variables, close all

% surface mesh of a sphere
Ps = trisphere(.4);

% tetrahedralization
DT = DelaunayTri([Ps; 0 0 0]);
P = DT.X;
T = DT.Triangulation;
M = 2*ones(size(T,1),1);

T = [T; T+size(P,1)];
P = [P; [P(:,1)+4 P(:,2:3)]];
M = [M; M*3/2];

% data structure
primal = createprimal3d(P,T,M);
plotscalar3d(primal,'mesh','y');

% electric data
r0 = 1;
V0 = 1;

% the problem has normal derivative as unknowns
% G*du = H*u

% imposed potential on the surface
fsurf = findsurf3d(primal);
B = barycenter(primal.node,primal.face_node(fsurf,:));
usurf = ones(length(fsurf),1);
usurf(B(:,1) > 1.5) = -1;

F = primal.face_node(fsurf,:);
idx = faceoutward3d(primal,fsurf,[]);
F(idx,:) = fliplr(F(idx,:));

% % % calculate G with function fkern
fkern1 = @(irow,jcol)bem_Gcoeff(irow,jcol,primal.node,F,B);
fkernACA = @(irow,jcol)bem_Gcoeff_num(irow,jcol,primal.node,F,B);
fkern2 = @(irow,jcol)bem_Hcoeff(irow,jcol,primal.node,F,B);
% 
% [G,H] = bem_matrix3d(primal);
% G2 = fkern1(1:length(fsurf),1:length(fsurf));
% H2 = fkern2(1:length(fsurf),1:length(fsurf));

tic
Gk = hmtx_cluster(B,'eta',1.5,'Nmin',50);
Hk = Gk;
% Hk = hmtx_cluster(B,'eta',2,'Nmin',50);
toc
% hmtx_plot(Gk);
eps = 1e-4;
tic; Gk = hmtx_fill(Gk,eps,fkern1); toc
tic
Hk = hmtx_fill(Hk,eps,fkern2);
Hk = hmtx_fixdiag(Hk);
toc

tic; Gkinv = hmtx_inv(Gk); toc
I = hmtx_mult(Gk,Gkinv);
If = hmtx_full(I);
figure,plot(diag(If.M))
Iz = If.M-diag(diag(If.M));
figure,plot(Iz(:))

% rhs
b = hmtx_HxM(Hk,usurf);

return

% solution
% dusurf = G\b;
errMAX = 1e-10;
% M = diag(sqrt(diag(G)));
% [dusurf,flag,relres,iter,resvec] = gmres(G,b,[],errMAX,500);
afun = @(x)Axfun(x,Gk);
[dusurf,flag,relres,iter,resvec] = gmres(afun,b,[],errMAX,500);
% 
% convergence
h = figure;
plot(1:iter(2)+1,log10(relres*resvec/resvec(end)),'r-');
hold on
plot(1:iter(2)+1,log10(errMAX)*ones(1,iter(2)+1),'b--'); % convergence tolerance
axis([0 iter(2)+1 log10(errMAX/10) 1])
xlabel('iteration #')
ylabel('log10(|| Ax-b ||/|| b ||)');
figuresetup(h);

% preallocate
u = zeros(primal.face_num,1);
du = zeros(primal.face_num,1);

% solution projected to volumic data structure
u(fsurf) = usurf;
du(fsurf) = dusurf;

% color maps of solution
plotscalar3d(primal,'field',u,'mesh','y','title','solution: u');
plotscalar3d(primal,'field',du,'mesh','y','title','solution: du');