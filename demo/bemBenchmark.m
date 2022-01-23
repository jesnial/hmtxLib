clear variables, close all

%% H-matrix params
% tolerance of the ACA method
tol = 1e-4;
% maximum rank
kMax = Inf;
% admissibility criterion
eta = 4;
% minimum points per cluster
Nmin = 32;      


%% geometric params
% radius of central sphere
rC = 3;
% number of satellites
nS = 10;
% radius of satellite sphere
rS = 1;
% radius of orbit
rO = 5;

% positive potential object code code
iObjPlus = 2;
% negative potential object code code
iObjMinus = 3;

% potential
V0 = 5;

% dimensions
hSize = .4;


%% geometry creation
% central
[P,T] = triSphere(hSize);
P = P*rC;
M = iObjPlus*ones(size(T,1),1);
% satellite
[PS,TS] = triSphere(hSize);
PS = PS*rS;

% load satellite
for i = 1:nS
    % angle
    alpha = 2*pi/nS*(i-1);
    % center
    C = [rO*cos(alpha) rO*sin(alpha) 0*alpha];
    % update
    T = [T; TS+size(P,1)];
    P = [P; bsxfun(@plus,rS*PS,C)];
    if mod(i,2) == 0
        M = [M; iObjPlus*ones(size(TS,1),1)];
    else
        M = [M; iObjMinus*ones(size(TS,1),1)];
    end
end

primal = createPrimal2d(P,T,M);
plotScalar2d(primal);
light

fprintf('Number of faces: %d\n',size(primal.Fac2Nod,1));

% target poits for the calulation of coefficients
Q = barycenter(P,T);
% [G,~] = bemMatrix2dm(P,T,Q);

% rhs
b = zeros(size(M));
b(M == iObjPlus) = +V0;
b(M == iObjPlus) = -V0;

% tic
% [G0,H0] = bemMatrix2d(primal);
% toc
% % G kernel function (from full matrix)
% fkern = @(irow,jcol)G0(irow,jcol);

% matrix
fkern = @(i,j)bemMatrix2dm(P,T(j,:),Q(i,:));
tic
G = hmtx_cluster(Q,'eta',eta,'Nmin',Nmin);
toc
tic
G = hmtx_fill(G,tol,fkern);
toc

% solve
tic
[L,U] = hmtx_lu(G);
toc
tic
y = hmtx_lsolve(L,b);
q = hmtx_usolve(U,y);
toc
plotScalar2d(primal,'field',q)

% solution: q = G\b;
errMAX = 1e-10;
iterMAX = 50;

AxFun = @(x)hmtx_HxM(G,x);
tic
[q,flag0,relres0,iter0,resvec0] = gmres(AxFun,b,[],errMAX,iterMAX,[]);
toc
