clear variables, close all;

% number of spheres
N = 3;

% params
IWire = 100;
rCoil = 20e-3;
nCoil = 50;
f = 85e3;

% sphere1
rSphere1 = 15e-3;
objSphere1 = 2;
hSizeSphere1 = 6e-3;
cSphere1 = [0 0 0];

% material
MuRSphere(1) = 1000;
MuRSphere(2) = 1000;
SigmaSphere(1) = 0;
SigmaSphere(2) = 1e6;

% core geometry and mesh
[P0,T0] = triSphere(hSizeSphere1/rSphere1);
P0 = bsxfun(@plus,P0*rSphere1,cSphere1);
M0 = ones(size(T0,1),1);

P = [];
T = [];
M = [];

for i = 1:N
    P = [P; P0(:,1:2) P0(:,3)+2.1*rSphere1*((i-1)-(N-1)/2)];
    T = [T; (i-1)*size(P0,1)+T0];
    M = [M; (i-1)+M0];
end

% create data structure
primal = createPrimal2d(P,T,M);
mat = setMaterialProperties(M);
mat(1).MuR = MuRSphere(1);
mat(2).MuR = MuRSphere(1);
mat(3).MuR = MuRSphere(2);
mat(1).Sigma = SigmaSphere(1);
mat(2).Sigma = SigmaSphere(1);
mat(3).Sigma = SigmaSphere(2);

clear P T M

h = plotScalar2d(primal,'mesh','y');

% coil
A = linspace(0,2*pi,nCoil).';
C0 = [rCoil*cos(A) rCoil*sin(A) 0*A];

PCoil1 = C0(1:end-1,:);
PCoil2 = C0(2:end,:);

plotLine3d(PCoil1,PCoil2,'color',[0 0 1],'handle',h);
print('-depsc2','-r600','object.eps');

% coil contribution
HWireFun = @(Q)segmentMagneticField3d(PCoil1,PCoil2,IWire,Q);

% solution
%tic; [phi,DphiDn] = bemSibcMagneticSolve2d(primal,mat,f,HWireFun,'method','constant','preconditioner','n'); toc
tic; [phi,DphiDn] = bemSibcMagneticSolve2d4(primal,mat,f,HWireFun,'method','constant','preconditioner','n'); toc
size(mat)
primal
% postprocessing
HStruct = @(Q)bemVectorField2d(primal,phi,DphiDn,Q,'verbose','n');
HFun = @(Q)HWireFun(Q)+HStruct(Q);

% field points
nQ = 1501;
rQ = (rCoil+rSphere1)/2;
Q = [rQ*ones(nQ,1) zeros(nQ,1) linspace(min(primal.Node(:,3)),max(primal.Node(:,3)),nQ).'];
figure(h), hold on
plot3(Q(:,1),Q(:,2),Q(:,3),'.');

H = HFun(Q);
H0 = HWireFun(Q);

s = sqrt(sum((Q-Q(1,:)).^2,2));
h = figure;
subplot(2,1,1); hold on
plot(s,real(H(:,1)));
plot(s,imag(H(:,1)));
subplot(2,1,2); hold on
plot(s,real(H(:,3)));
plot(s,imag(H(:,3)));
print('-depsc2','-r600','res1.eps');

% load femm
Ht = load('bem3spheres85kHzHt.txt');
Hn = load('bem3spheres85kHzHn.txt');

subplot(2,1,1);
plot(Hn(:,1),-Hn(:,3:end),':');
l = legend('$\mathcal{H}$-matrix real','$\mathcal{H}$-matrix imag','FEMM real', 'FEMM imag');
%l = legend('Sibc-BEM real', 'Sibc-BEM imag', 'FEMM real', 'FEMM imag');
set(l, 'interpreter', 'latex');
xlabel('Radial distance (m)') 
ylabel('Magnetic flux density (mT)') 
subplot(2,1,2);
plot(Ht(:,1),Ht(:,3:end),':');
l = legend('$\mathcal{H}$-matrix real','$\mathcal{H}$-matrix imag','FEMM real', 'FEMM imag');
%l = legend('Sibc-BEM real', 'Sibc-BEM imag', 'FEMM real', 'FEMM imag');
set(l, 'interpreter', 'latex');
xlabel('Radial distance (m)') 
ylabel('Magnetic flux density (mT)')
print('-depsc2','-r600','res2.eps');


