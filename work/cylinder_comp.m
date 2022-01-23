%% cylinder
% BEM formulation of two cylindrical electrodes with imposed potentials 'u'
% the problem has normal derivative of potential 'du' as unknowns
% G*du = H*u = b

clear variables, close all

%% parameters
r1 = .5;     % inner radius
r2 = .7;     % outer radius
dz = 2;      % heigth of cylinder
V0 = 2;      % voltage applied between the eletrodes
eps = 1e-4;  % tolerance of the ACA method
tol = 1e-10; % tolerance for GMRES
% h = 0.035;    % mesh size;

h = .1:-.005:.02;  % mesh size;

nh = length(h);
nF = zeros(nh,1);
tassembly1 = zeros(nh,1);
MB1 = zeros(nh,1);
tproduct1 = zeros(nh,1);
tassembly2 = zeros(nh,1);
MB2 = zeros(nh,1);
tproduct2 = zeros(nh,1);

for i = 1:nh
    %% data structure: calculation of nodes P, faces F, barycenters B
    % number of divisions
    na1 = round(2*pi*r1/h(i));
    na2 = round(2*pi*r2/h(i));
    nz = round(dz/h(i));
    
    % surface mesh of inner cylinder
    [P1,F1] = thincylinder(r1,-dz/2,dz/2,na1,nz);
    
    % surface mesh of outer cylinder
    [P2,F2] = thincylinder(r2,-dz/2,dz/2,na2,nz);
    
    % concatenate structure
    F = [F1; F2+size(P1,1)];
    P = [P1; P2];
    
    % barycenters of triangles
    B = barycenter(P,F);
    
    % random vector
    b = rand(size(B,1),1);
    
    nF(i) = size(F,1);
    
    %% H matrix
    % calculate G with function fkern
    fkern1 = @(irow,jcol)bem_Gcoeff(irow,jcol,P,F,B);
    
    % matrix population
    t0 = tic;
    Gk = hmtx_cluster(B,'eta',2,'Nmin',50);
    Gk = hmtx_fill(Gk,eps,fkern1);
    tassembly1(i) = toc(t0);
    
    % matrix memory occupation
    MB1(i) = hmtx_Mbyte(Gk);
    
    % matrix-vector product
    t0 = tic;
    for j = 1:10
        x1 = hmtx_HxM(Gk,b);
    end
    tproduct1(i) = toc(t0);
    
    clear Gk
    
    %% full matrix
    % matrix population
    t0 = tic;
    G = bem_Gcoeff(1:size(F,1),1:size(F,1),P,F,B);
    tassembly2(i) = toc(t0);
    
    % matrix memory occupation
    MB2(i) = getfield(whos('G'),'bytes')/2^20;
    
    % matrix-vector product
    t0 = tic;
    for j = 1:10
        x2 = G*b;
    end
    tproduct2(i) = toc(t0);
    clear G
end

save('scalability.mat','nF','tassembly1','tproduct1','MB1','tassembly2','tproduct2','MB2')
% % tassembly
% figure; hold on
% id1 = tassembly1 ~= 0;
% plot(nF(id1),tassembly1(id1),'ro');
% id2 = tassembly2 ~= 0;
% plot(nF(id2),tassembly2(id2),'bo');
% legend('H-matrix','full-matrix')
% xlabel('number of elements'); ylabel('computational time (s)')
% title('matrix assembly')
% 
% % tproduct
% figure; hold on
% id1 = tproduct1 ~= 0;
% plot(nF(id1),tproduct1(id1),'ro');
% id2 = tproduct2 ~= 0;
% plot(nF(id2),tproduct2(id2),'bo');
% legend('H-matrix','full-matrix')
% xlabel('number of elements'); ylabel('computational time (s)')
% title('matrix vector product (10x)')
% 
% % memory
% figure; hold on
% id1 = MB1 ~= 0;
% plot(nF(id1),MB1(id1),'ro');
% id2 = MB2 ~= 0;
% plot(nF(id2),MB2(id2),'bo');
% legend('H-matrix','full-matrix')
% xlabel('number of elements'); ylabel('memory (MB)')
% title('memory allocation')