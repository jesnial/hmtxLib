% SVD: reconstruction of potential

clear variables, close all

% number of sources
N = 1000;
% number of targets
M = 2000;

% radius of sources
R0 = 1;
% inner radius of targets
R1 = 0;
% outer radius of targets
R2 = 10;

% desired tolerance
eps = logspace(-1,-10,10);

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

% plot
h1 = figure; hold on, axis equal, view([1 1 1])
plot3(P(:,1),P(:,2),P(:,3),'o','MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor','b');
plot3(Q(:,1),Q(:,2),Q(:,3),'o','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r');

% calculate 1/r matrix (assume point charge q = 1)
D = distance(Q,P);
% kernel matrix K
K = 1./D;

% exact potential
Vex = K*q;

% singular value decomposition
[U,S,V] = svd(K);
sigma = diag(S);

% get truncation
% sigma in increasing order
sigma2 = sort(sigma);
normrKr2 = cumsum(sigma2.^2);

% norm of K
normK2 = norm(K,'fro')^2;
h2 = figure; hold on
%plot(1:M,Vex,'r')

for i = 1:length(eps)
    ir = find(normrKr2 <= eps(i)^2*normK2,1,'last');
    r = length(sigma)-ir+1;
    
    nr(i) = r;
    Kr = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    Vr(:,i) = Kr*q;
end
plot(sqrt(sum(Q.^2,2)),log10(abs(Vr-repmat(Vex,1,length(eps)))./repmat(Vex,1,length(eps))))
xlabel('curvilinear coordinate');
ylabel('log10 of relative error on potential');
legendCell = cellstr(num2str([eps' nr'], 'eps = %5.0e, r = %-d'));
legend(legendCell);
figuresetup(h2);
