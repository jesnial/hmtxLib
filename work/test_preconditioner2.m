clear variables, close all

% electric data
r0 = 1;
V0 = 1;

% Hmatrix params
eps = 1e-6;
eta = 3;
Nmin = 32;

% gmres params
errMAX = 1e-10;
iterMAX = 40;

% surface mesh of a sphere
[P,F] = triSphere(.2);

% normalize radius to r0
r = sqrt(sum(P.^2,2));
P = r0*P./repmat(r,1,3);

figure; patch('Faces',F,'Vertices',P,'FaceColor','r','EdgeColor','k');
view([1 1 1]); axis equal; colorbar

% the problem has normal derivative as unknowns
% G*du = H*u

% % calculate G with function fkern
B = barycenter(P,F);
Gkern = @(irow,jcol)bem_Gcoeff(irow,jcol,P,F,B);


tic; G = hmtx_cluster(B,'eta',eta,'Nmin',Nmin); toc
H = hmtx_copystruct(G);

tic; G = hmtx_fill(G,eps,Gkern); toc

% rhs
b = V0*ones(G.nrow,1);


afun = @(x)hmtx_mvm(G,x,zeros(G.nrow,1),1,1);
tic; D = hmtx_blkdiagprecond(G); toc
tic; [L,U] = hmtx_ssor(G); toc

Mfun = @(x)hmtx_mvm(D,x,zeros(G.nrow,1),1,1);
Mfun2 = @(x)hmtx_usolve(U,hmtx_lsolve(L,x));
M1fun = @(x)hmtx_lsolve(L,x);
M2fun = @(x)hmtx_usolve(U,x);

tic; [x1,flag1,relres1,iter1,resvec1] = gmres(afun,b,[],errMAX,iterMAX); toc
tic; [x2,flag2,relres2,iter2,resvec2] = gmres(afun,b,[],errMAX,iterMAX,Mfun); toc 
tic; [x3,flag3,relres3,iter3,resvec3] = gmres(afun,b,[],errMAX,iterMAX,M1fun,M2fun); toc
tic; [x4,flag4,relres4,iter4,resvec4] = gmres(afun,b,[],errMAX,iterMAX,Mfun2); toc 

% convergence
h = figure;
hold on
plot(1:length(resvec1),log10(relres1*resvec1/resvec1(end)),'r-');
plot(1:length(resvec2),log10(relres2*resvec2/resvec2(end)),'g-');
plot(1:length(resvec3),log10(relres3*resvec3/resvec3(end)),'b-');
plot(1:length(resvec4),log10(relres4*resvec4/resvec4(end)),'c-');
axis([0 iterMAX log10(errMAX/10) 1])
xlabel('iteration #')
ylabel('log10(|| Ax-b ||/|| b ||)');
legend('no precond','block diagonal','ssor');
plot(1:iterMAX,log10(errMAX)*ones(1,iterMAX),'k--'); % convergence tolerance

figuresetup(h);

% color maps of solution
h2 = figure; patch('Faces',F,'Vertices',P,'FaceVertexCData',c,'FaceColor','flat','EdgeColor','k');
view([1 1 1]); axis equal; colorbar
figuresetup(h2);