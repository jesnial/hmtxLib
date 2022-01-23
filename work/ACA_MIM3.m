clear variables, close all

%% parameters
sigma = 59.5e6;  % electrical conductivity
delta = 2e-3;    % thickness
f = 1000;        % frequancy

eps = 1e-4;  % tolerance of the ACA method
eta = .5;    % threshold for admissibility
Nmin = 32;   % minimum number of points in a leaf of the cluster tree
tol = 1e-8;  % tolerance for GMRES

%% mesh data structure
load('xbencha.mat','I','Ctrack');
I0 = I;
[P,T,M,track] = createtrack(Ctrack,I0,0);
[P,T,M] = refinesurf(P,T,M);
[P,T,M] = refinesurf(P,T,M);

%% create data structure
primal = createprimal2d(P,T,M);
plotscalar2d(primal,'mesh','y');
mat = setmatproperties(M);
clear P T M Ctrack I

%% material properties
mat.sigma(1:mat.nmat) = sigma;
mat.delta(1:mat.nmat) = delta;

%% define boundary
[idfix,xfix,idfloat] = boundaryconditions(primal,track);

% %% calculate coefficients
% % calculate inductance matrix
% % tic; L = newinductancematrixf(primal.node,int32(primal.face_node),int32(1:primal.node_num)',int32(1:primal.node_num)'); toc
% L = inductancematrix(primal,mat);
% % calculate resistance
% R = resistancematrix(primal,mat);
% 
% %% assemble matrix and rhs
% A = 1i*2*pi*f*L;
% b = zeros(size(A,1),1);
% 
% %% solution
% x1 = problemsolve(A,b,idfix,xfix,idfloat);
% % plot stream function
% plotscalar2d(primal,'field',real(x1),'mesh','y','isolevel',20,'expr','x<0 & y <0');
% 
% %% postprocessing
% Jt = postJ(primal,mat,x1);
% Ji = imag(Jt);
% Jr = real(Jt);
% plotvect2d(primal,Jr);
% 
% %% power
% PW1 = sum(powerlosses(primal,mat,Jt))


%% ACA-Hmatrix approach
clear mim_Lcoeff mim_Lcoeff_num  % necessary to erease the persistent variable inside the funciton
fkern1 = @(irow,jcol)mim_Lcoeff(irow,jcol,primal.node,primal.face_node);
fkern2 = @(irow,jcol)mim_Lcoeff_num(irow,jcol,primal.node,primal.face_node);

tic; Lk = hmtx_cluster(primal.node,'eta',eta,'Nmin',Nmin); toc

% partition
h = plotscalar2d(primal,'mesh','y');
c = hmtx_getcluster(Lk);
figure(1); hold on
for i = 1:length(c)
    col = rand(1,3);
    plot3(primal.node(c{i},1),primal.node(c{i},2),primal.node(c{i},3),'o',...
            'MarkerEdgeColor',col,'MarkerFaceColor',col,'MarkerSize',10);
end

% filling matrix
tic; Lk = hmtx_fill(Lk,eps,fkern1,fkern2); toc
b0 = zeros(Lk.nrow,1);
x2 = problemsolve_Hmatrix(Lk,b0,f,idfix,xfix,idfloat,'errMAX',tol);

plotscalar2d(primal,'field',real(x2),'mesh','y','isolevel',20,'expr','x<0 & y <0');

%% postprocessing
Jt = postJ(primal,mat,x2);
Ji = imag(Jt);
Jr = real(Jt);
plotvect2d(primal,Jr,'obj',1,'scale',1e-6);

%% power
PW2 = sum(powerlosses(primal,mat,Jt))


% figure(1); hold on
% z = ones(primal.face_num,1);
% for i = 1:length(c)
%     icell = find(sum(ismember(primal.face_node,c{i}),2) == 3);
%     z(icell) = i+1;
% end
% primal.obj = z;
% plotscalar2d(primal,'mesh','y','obj',23);
