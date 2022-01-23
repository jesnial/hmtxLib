% SVD: plot of singular values at different distances between sources and targets 

clear variables, close all

% number of sources
N = 1000;
% number of targets
M = 2000;

% radius of sources
R0 = 1;
% inner radius of targets
Rin = [0 1 2 9];
% outer radius of targets
Rout = [1 2 3 10];

% desired tolerance
eps = 1e-4;

% source points
r = R0*rand(N,1);
theta = pi*rand(N,1);
phi = 2*pi*rand(N,1);
P = [r.*sin(theta).*cos(phi) r.*sin(theta).*sin(phi) r.*cos(theta)];

for i = 1:length(Rin)
    
    % target points
    r = Rin(i)+(Rout(i)-Rin(i))*rand(M,1);
    theta = pi*rand(M,1);
    phi = 2*pi*rand(M,1);
    Q = [r.*sin(theta).*cos(phi) r.*sin(theta).*sin(phi) r.*cos(theta)];
    
    % plot
    h1 = figure; hold on, axis equal, view([1 1 1])
    plot3(P(:,1),P(:,2),P(:,3),'o','MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor','b');
    plot3(Q(:,1),Q(:,2),Q(:,3),'o','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r');
    figuresetup(h1); box off
    print('-dpdf',strcat('sphere',num2str(Rin(i)))); 
    
    % calculate 1/r matrix (assume point charge q = 1)
    D = distance(P,Q);
    % kernel matrix K
    K = 1./D;
    
    % perform SVD
    [U,S,V] = svd(K);
    sigma(:,i) = diag(S);
end

h2 = figure;
plot(log10(sigma));
xlabel('singular value No.');
ylabel('log10 singular value');
legend('Rin = 0, Rout = 1','Rin = 1, Rout = 2','Rin = 2, Rout = 3','Rin = 9, Rout = 10');
figuresetup(h2);
print('-dpdf',strcat('svd',num2str(Rin(i)))); 

