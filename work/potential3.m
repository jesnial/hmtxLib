% ACA with full pivoting: reconstruction of potential 

clear variables, close all

% number of sources
N = 4;
% number of targets
M = 5;

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
tic
Vex = K*q;
toc
% norm of K
normK = norm(K,'fro')^2;
h2 = figure; hold on
%plot(1:M,Vex,'r')

for i = 1:length(eps)
    % initialization
    R = K;
    err = 1;
    
    U = [];
    V = [];
    while err > eps(i)
        % find pivot
        [maxK,imaxR] = max(R(:));
        [imax,jmax] = ind2sub(size(K),imaxR);
        % normalizing constant
        gamma = 1./R(imaxR);
        % new functions
        u = gamma*R(:,jmax);
        v = R(imax,:);
        % new residual
        R = R-u*v;
        % new approximation
        U = [U u];
        V = [V; v];
        normR = norm(R,'fro');
        err = normR/normK;
    end
    
    nr(i) = size(U,2);
    tic
    Vr(:,i) = U*(V*q);
    toc
end
plot(sqrt(sum(Q.^2,2)),log10(abs(Vr-repmat(Vex,1,length(eps)))./repmat(Vex,1,length(eps))))
xlabel('curvilinear coordinate');
ylabel('log10 of relative error on potential');
legendCell = cellstr(num2str([eps' nr'], 'eps = %5.0e, r = %-d'));
legend(legendCell);
figuresetup(h2);
