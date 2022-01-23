function [idx1,idx2] = bisect3d(P)

% BISECT3D partitions the points P into two disjoint sets
%
% USE:
% [idx1,idx2] = bisect3d(P)
%
% INPUTS:
% 'P': point coordinates
%
% OUTPUTS:
% 'idx1': indices of first cluster
% 'idx2': indices of second cluster
%
% NOTE:
% See S. Borm, L. Grasedyck, W. Hackbusch, "THierarchical matrices", 2005,
% pp. 26--29
%
% VERSION:
% Date: 08.01.2013
% Copyright(C) 2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:

% find coordinate with maximum extent
Pmin = min(P,[],1);
Pmax = max(P,[],1);
delta = Pmax-Pmin;
[~,imax] = max(delta);

% midpoint
Pmid = .5*(Pmax(imax)+Pmin(imax));

% clustering
idx1(:,1) = find(P(:,imax) > Pmid);
idx2(:,1) = find(P(:,imax) <= Pmid);
