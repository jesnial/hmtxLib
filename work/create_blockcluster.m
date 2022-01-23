clear variables, close all

Ps = trisphere(.08);
% Ps(:,3) = Ps(:,3)*2;

% tetrahedralization
DT = DelaunayTri([Ps; 0 0 0]);
P = DT.X;
T = DT.Triangulation;

% T = [T; T+size(P,1)];
% P = [P; P(:,1)+4 P(:,2:3)];
M = 2*ones(size(T,1),1);

% data structure
primal = createprimal3d(P,T,M);
plotscalar3d(primal,'mesh','y');

% surface faces
fsurf = findsurf3d(primal);

% barycenters
P = primal.node;
T = primal.face_node(fsurf,:);
B = barycenter(P,T);

eps = 1e-4;
fkern = @(irow,jcol)one_dividedby_r(irow,jcol,B,B);

tic
A = hmtx_cluster(B,'eta',.1,'Nmin',100);
toc
tic
A = hmtx_fill(A,fkern,eps);
toc
% h = hmtx_plot(A);

% tic
% Af = fkern(1:A.nrow,1:A.ncol);
% toc
% 
b = rand(length(fsurf),1);
% tic; xf = Af*b; toc;
tic; x = hmtx_Ax(A,b); toc;
figure,plot(xf,'ro');hold on, plot(x,'b*')
figure,plot(log10(abs(x-xf)./xf))
% h = plothmatrix(A);

% for i = 1:nlev
%     % create all possible cluster pairs
%     pairs = nchoosek(1:2^i,2);