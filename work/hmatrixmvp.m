%% hmtxmvp test complexity of matrix vector product

clear variables, close all

%% data structure: calculation of nodes P, faces F, barycenters B

% icosahedron
phi = (1+sqrt(5))/2;
P = [
    0 0 0
    0  1  phi
    0 -1  phi
    0  1 -phi
    0 -1 -phi
    1  phi 0
    -1  phi 0
    1 -phi 0
    -1 -phi 0
    phi 0  1
    phi 0 -1
    -phi 0  1
    -phi 0 -1
    ];

% normalize radius to 1
P = P./sqrt(1+phi^2);

% triangulate
DT = DelaunayTri(P);

% surface boundary
FB = freeBoundary(DT);

nstep = 7;

% preallocation
tclus = zeros(nstep,1);
tfill = zeros(nstep,1);
tcomp = zeros(nstep,1);
tmvp1 = zeros(nstep,1);
tmvp2 = zeros(nstep,1);
nunk = zeros(nstep,1);
stor1 = zeros(nstep,1);
stor2 = zeros(nstep,1);

for i = 1:nstep
    
    % get face barycenter
    B = barycenter(DT.X,FB);
    
    nunk(i) = size(B,1);
    
    % get points and faces
    P = DT.X(2:end,:);
    F = FB-1;
    
    % color maps of solution
    figure; patch('Faces',F,'Vertices',P,'FaceColor','r','EdgeColor','k');
    view([1 1 1]); axis equal; colorbar
    
    %% matrix calculation
    % calculate single layer BEM matrix with function fkern
    fkern1 = @(irow,jcol)bem_Gcoeff(irow,jcol,P,F,B);
    
    t0 = tic;
    Gk = hmtx_cluster(B,'eta',1.5,'Nmin',32);
    tclus(i) = toc(t0);
    
    stor1(i) = hmtx_Mbyte(Gk);

    % populate Gk
    t0 = tic;
    Gk = hmtx_fill(Gk,eps,fkern1);
    tfill(i) = toc(t0);
% 
%     % vector
%     x = ones(Gk.nrow,1);
% 
%     % rhs
%     t0 = tic;
%     for j = 1:10
%         b = hmtx_HxM(Gk,x);
%     end
%     tmvp1(i) = toc(t0);
% 
%     t0 = tic;
%     Gk = hmtx_compress(Gk);
%     tcomp(i) = toc(t0);
%     
%     stor2(i) = getfield(whos('Gk'),'bytes');
%     
%     
%     % rhs
%     t0 = tic;
%     for j = 1:100
%         b = hmtx_HxM(Gk,x);
%     end
%     tmvp2(i) = toc(t0);
    
    % update icosahedron
    % normalize to unit sphere
    r = sqrt(sum(B.^2,2));
    B = B./repmat(r,1,3);
    
    % triangulate
    DT = DelaunayTri([DT.X; B]);
    
    % surface boundary
    FB = freeBoundary(DT);
    
end