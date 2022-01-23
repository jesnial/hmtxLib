clear variables%, close all

Ps = triSphere(.2);
% Ps(:,3) = Ps(:,3)*2;

% tetrahedralization
DT = DelaunayTri([Ps; 0 0 0]);
P = DT.X;
T = DT.Triangulation;

T = [T; T+size(P,1)];
P = [P; P(:,1)+4 P(:,2:3)];
M = 2*ones(size(T,1),1);

% convert to surface
[P,T,M] = volume2surface(P,T,M);

B = barycenter(P,T);

v1 = P(T(:,2),:)-P(T(:,1),:);
v2 = P(T(:,3),:)-P(T(:,1),:);
A = .5*cross(v1,v2,2);
A = sqrt(sum(A.^2,2));

nlev = 7;

% default value
cluster = zeros(size(B,1),nlev);

for i = 1:nlev
    % cluster codes
    iclus = unique(cluster(:,1:i),'rows');
    for j = 1:length(iclus)
        idx = find(ismember(cluster(:,1:i),iclus(j,:),'rows'));
        tic; [idx1,idx2] = bisect3d(B(idx,:)); toc
        col = zeros(size(idx));
        col(idx2) = 1;
        cluster(idx,i) = col;
        
%         % get unique code for different clusters
%         [~,code] = ismember(cluster,unique(cluster,'rows'),'rows');
%         plotscalar3d(primal,'mesh','y','field',code);

    end
end

% get unique code for different clusters
[~,code] = ismember(cluster,unique(cluster,'rows'),'rows');

pp = createPrimal2d(P,T,code+1);
plotScalar2d(pp,'mesh','y');