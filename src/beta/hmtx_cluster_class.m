function A = hmtx_cluster_class(P,varargin)

% HMTX_CLUSTER creates the supermatrix structure of cluster of points
%
% USE:
% [A,idx] = hmtx_cluster(P,varargin)
%
% INPUTS:
% 'P': cluster of points
% 'varargin': optional inputs
%   * 'eta': parameter for the admissibility
%   * 'Nmin': minimum number of points in a leaf of the clustertree
%
% OUTPUTS:
% 'A': supermatrix structure
%
% NOTE:
% Default Nmin = 32: see S. Borm, L. Grasedyck, W. Hackbusch, "Introduction
% to Hierachical Matrices with Application", 2002, pp. 4
%
% VERSION:
% Date: 09.01.2013
% Copyright(C) 2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:

% set defaults
eta = .1;
Nmin = 32;

% check varargin
for i = 1:2:length(varargin)
    vin = varargin{i};
    if strcmpi(vin,'eta')
        eta = varargin{i+1};
    elseif strcmpi(vin,'Nmin')
        Nmin = varargin{i+1};
    end
end

% create empty H-matrix
A = Hmatrix('supermatrix',1:size(P,1),1:size(P,1));

% create matrix cluster structure
A = matrixcluster(A,P,(1:A.nrow)',(1:A.ncol)',eta,Nmin);

end


% function called recursively
function A = matrixcluster(A,P,idx1,idx2,eta,Nmin)

% 31.12.2012

% 03.01.2013: added index values
% 10.01.2012: changed call to CLUSTERADMISSIBILITY
% 17.01.2013: uses class Hmatrix

% here only 'supermatrix' type can enter

% bisection
[idxM{1},idxM{2}] = bisect3d(P(idx1,:));
[idxN{1},idxN{2}] = bisect3d(P(idx2,:));

% local-to-global
idxM{1} = idx1(idxM{1});
idxM{2} = idx1(idxM{2});
idxN{1} = idx2(idxN{1});
idxN{2} = idx2(idxN{2});

% create structure
for im = 1:2
    for jn = 1:2
        
        % check admissibility of blocks
        if clusteradmissibility(P(idxM{im},:),P(idxN{jn},:),eta);
            A.M{im,jn} = Hmatrix('rkmatrix',idxM{im},idxN{jn});
        elseif length(idxM{im}) < Nmin || length(idxN{jn}) < Nmin
            A.M{im,jn} = Hmatrix('fullmatrix',idxM{im},idxN{jn});
        else
            A.M{im,jn} = Hmatrix('supermatrix',idxM{im},idxN{jn});
            A.M{im,jn} = matrixcluster(A.M{im,jn},P,A.M{im,jn}.irow,A.M{im,jn}.jcol,eta,Nmin);
        end
        
    end
end

% return reordered indices
A.irow = [A.M{1,1}.irow; A.M{2,1}.irow];
A.jcol = [A.M{1,1}.jcol; A.M{1,2}.jcol];

end

