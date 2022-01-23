function out = clusteradmissibility(P1,P2,eta)

% CLUSTERADMISSIBILITY checks if two cluster are admissible
%
% USE:
% out = clusteradmissibility(P1,P2,eta)
%
% INPUTS:
% 'P1': first cluster of points
% 'P2': second cluster of points
% 'eta': parameter for the admissibility
%
% OUTPUTS:
% 'out': logical variable (1: admissible, 0: not admissible)
%
% NOTE:
%
% VERSION:
% Date: 29.12.2012
% Copyright(C) 2012-2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 03.01.2013: added exact admissibility check
% 03.01.2013: added bounding-box check
% 10.01.2013: use of two different clusters as input

method = 3;

if method == 1
    % diameters of clusters
    D1 = max(reshape(distance(P1,P1),size(P1,1)*size(P1,1),1));
    D2 = max(reshape(distance(P2,P2),size(P2,1)*size(P2,1),1));
    % distance between clusters
    d = min(reshape(distance(P1,P2),size(P1,1)*size(P2,1),1));
elseif method == 2
    % centers of clusters
    C1 = sum(P1,1)/size(P1,1);
    C2 = sum(P2,1)/size(P2,1);
    % diameters of clusters
    D1 = 2*max(sqrt(sum(bsxfun(@minus,P1,C1).^2,2)));
    D2 = 2*max(sqrt(sum(bsxfun(@minus,P2,C2).^2,2)));
    % distance between clusters
elseif method == 3
    % cluster bounding boxes
    minP1 = min(P1,[],1);
    maxP1 = max(P1,[],1);
    minP2 = min(P2,[],1);
    maxP2 = max(P2,[],1);
    
    % origin
    Z = zeros(size(minP1));
    % diameters of clusters
    D1 = sqrt(sum((maxP1-minP1).^2));
    D2 = sqrt(sum((maxP2-minP2).^2));
    % distance between clusters
    d = sqrt(sum(max(Z,minP1-maxP2).^2+max(Z,minP2-maxP1).^2,2));
end
out = min([D1 D2]) < eta*d;

end

