% ACA with partial pivoting: reconstruction of potential

clear variables%, close all

% number of sources
N = 10000;
% number of targets
M = 20000;

% radius of sources
R0 = 1;
% inner radius of targets
R1 = 0;
% outer radius of targets
R2 = 10;

% desired tolerance
eps = 1e-8; %logspace(-1,-10,10);

% source points
r = R0*rand(N,1);
theta = pi*rand(N,1);
phi = 2*pi*rand(N,1);
P = [r.*sin(theta).*cos(phi) r.*sin(theta).*sin(phi) r.*cos(theta)];

% charges
q = ones(N,1);

% target points
s = sqrt(3)*linspace(R1,max(R2),M)';
Q = [s s s];

% % plot
% h1 = figure; hold on, axis equal, view([1 1 1])
% plot3(P(:,1),P(:,2),P(:,3),'o','MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor','b');
% plot3(Q(:,1),Q(:,2),Q(:,3),'o','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r');

% kernel matrix K
tic
K = 1./distance(Q,P);
toc

% exact potential
tic
Vex = K*q;
toc
% norm of K
normK = norm(K,'fro')^2;
%plot(1:M,Vex,'r')

nrmax = max(round(min([M,N])/10),100);
for i = 1:length(eps)

    krow = @(irow)kernelrow(irow,P,Q);
    kcol = @(jcol)kernelcol(jcol,P,Q);
    tic
    [U,VT] = aca_old(krow,kcol,M,N,eps(i));
    toc
    
    nr(i) = size(U,2);
    tic
    Vr(:,i) = U*(VT*q);
    toc
end
h2 = figure; hold on
plot(sqrt(sum(Q.^2,2)),log10(abs(Vr-repmat(Vex,1,length(eps)))./repmat(Vex,1,length(eps))))
xlabel('curvilinear coordinate');
ylabel('log10 of relative error on potential');
legendCell = cellstr(num2str([eps' nr'], 'eps = %5.0e, r = %-d'));
legend(legendCell);
figuresetup(h2);
