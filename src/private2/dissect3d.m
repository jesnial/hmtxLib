function [idx1,idx2] = dissect3d(P)

% DISSECT3D partition the points P into two disjoint sets
%
% USE:
% col = dissect3d(P)
%
% INPUTS:
% 'P': point coordinates
%    
% OUTPUTS:
% 'idx1': indices of first cluster
% 'idx2': indices of second cluster
%
% NOTE:
% See S. Rjasanow, O. Steinbach, "The fast solution of Boundary Integral
% Equations", Springer, 2007, pp. 109
%
% VERSION:
% Date: 17.12.2012
% Copyright(C) 2012-2020: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 18.12.2012: output: 0/1
% 18.12.2012: added check for max eigenvalue
% 29.12.2012: weights are optional
% 17.04.2013: bug fix in covariance matrix calculation
% 03.01.2020: uses built-in COV function

% centre of the cluster
Xc = mean(P,1);

% covariance matrix of the cluster
C = cov(P);

% eigenvalues and eigenvectors
[V,D] = eig(C);

% get max eigenvalue and corresponding eigenvector
lambda = diag(D);
[~,iMax] = max(abs(lambda));
vMax = V(:,iMax);

% separation
dist = bsxfun(@minus,P,Xc)*vMax;  % dot product

% clustering
idx1 = find(dist > 0);
idx2 = find(dist <= 0);

% % split
% dist = sort(dist);
% iDiv = round(length(dist)/2);
% idx1 = 1:iDiv;
% idx2 = iDiv+1:length(dist);

end

% if nargin == 1
%     g = ones(size(P,1),1);
% end
% 
% % mass of the cluster
% G = sum(g);
% 
% % centre of the cluster
% X = sum(bsxfun(@times,g,P),1)/G;
% 
% % covariance matrix of the cluster
% Xt = -bsxfun(@minus,P,X);
% C = cov(Xt);
% C = [
%     sum(g.*Xt(:,1).^2) sum(g.*Xt(:,1).*Xt(:,2)) sum(g.*Xt(:,1).*Xt(:,3))
%     sum(g.*Xt(:,2).*Xt(:,1)) sum(g.*Xt(:,2).^2) sum(g.*Xt(:,2).*Xt(:,3))
%     sum(g.*Xt(:,3).*Xt(:,1)) sum(g.*Xt(:,3).*Xt(:,2)) sum(g.*Xt(:,3).^2)
%     ];
% 
% % eigenvalues and eigenvectors
% [V,D] = eig(C);
% 
% % get max eigenvalue and corresponding eigenvector
% lambda = diag(D);
% [~,iMax] = max(abs(lambda));
% vMax = V(:,iMax);
% 
% % separation
% dist = Xt*vMax;  % dot product
% 
% % clustering
% idx1 = find(dist > 0);
% idx2 = find(dist <= 0);
% 
% % % split
% % dist = sort(dist);
% % iDiv = round(length(dist)/2);
% % idx1 = 1:iDiv;
% % idx2 = iDiv+1:length(dist);
% 
% end
