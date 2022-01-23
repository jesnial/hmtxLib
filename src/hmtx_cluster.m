function A = hmtx_cluster(P1,varargin)

% HMTX_CLUSTER creates the supermatrix structure of cluster of points
%
% USE:
% A = hmtx_cluster(P1,varargin)
%
% INPUTS:
% 'P1': cluster of points
% 'varargin': optional inputs
%   * 'eta': parameter for the admissibility (default, eta = 2)
%   * 'Nmin': minimum number of points in a leaf of the cluster tree (default, Nmin = 32)
%   * 'P2': second cluster of points (default, P2 = P1)
%   * 'kMax' : maximum rank (default: Inf: unconstrained rank)
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
% Copyright(C) 2013-2020: Fabio Freschi (fabio.freschi@polito.it)
%                         
% HISTORY:
% 17.01.2013: generalized for two clusters of points
% 23.01.2013: H-matrix is not initialized before calling MATRIXCLUSTER
% 03.01.2020: added warning for wrong optional inputs


% set defaults
eta = 2;
Nmin = 32;
P2 = P1;
iStart = 0;
jStart = 0;

% check varargin
for i = 1:2:length(varargin)
    vin = varargin{i};
    if strcmpi(vin,'eta')
        eta = varargin{i+1};
    elseif strcmpi(vin,'Nmin')
        Nmin = varargin{i+1};
    elseif strcmpi(vin,'P2')
        P2 = varargin{i+1};
    elseif strcmpi(vin,'iStart')
        iStart = varargin{i+1};
    elseif strcmpi(vin,'jStart')
        jStart = varargin{i+1};
    else
        warning('MATLAB:hmtx_cluster','Wrong VARARGIN parameter, skip ''%s = %s''\n',vin,num2str(varargin{i+1}));
    end
end

% create matrix cluster structure
A = matrixcluster(P1,(1:size(P1,1)).',P2,(1:size(P2,1)).',eta,Nmin,iStart,jStart);

end

function A = matrixcluster(P1,idx1,P2,idx2,eta,Nmin,iStart,jStart)

% MATRIXCLUSTER recursively creates the supermatrix structure starting from
% two clusters of points
%
% USE:
% A = matrixcluster(A,P1,idx1,P2,idx2,eta,Nmin)
%
% INPUTS:
% 'A': supermatrix structure
% 'P1': first cluster of points
% 'idx1': subset of P1 to be processed
% 'P2': second cluster of points
% 'idx2': subset of P2 to be processed
% 'eta': parameter for the admissibility
% 'Nmin': minimum number of points in a leaf of the clustertree
% 'kMax' : maximum rank (default: Inf: unconstrained rank)
%
% OUTPUTS:
% 'A': updated supermatrix structure
%
% NOTE:
% Function called recursively
% 
% VERSION:
% Date: 31.12.2012
% Copyright(C) 2012-2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 03.01.2013: added index values
% 10.01.2013: changed call to CLUSTERADMISSIBILITY
% 18.01.2013: refurbished routine
% 23.01.2013: uses HMTX_CREATE to create H-matrix structure

% check admissibility of blocks
if clusteradmissibility(P1(idx1,:),P2(idx2,:),eta)
    A = hmtx_create('rkmatrix',idx1+iStart,idx2+jStart);
elseif length(idx1) < Nmin || length(idx2) < Nmin
    A = hmtx_create('fullmatrix',idx1+iStart,idx2+jStart);
else
    A = hmtx_create('supermatrix',idx1+iStart,idx2+jStart);
    
    % dissection
    [idx1A,idx1B] = dissect3d(P1(idx1,:));
    [idx2A,idx2B] = dissect3d(P2(idx2,:));
    
    % local-to-global
    idx1A = idx1(idx1A);
    idx1B = idx1(idx1B);
    idx2A = idx2(idx2A);
    idx2B = idx2(idx2B);
    
    % recursive call
    A.M{1,1} = matrixcluster(P1,idx1A,P2,idx2A,eta,Nmin,iStart,jStart);
    A.M{1,2} = matrixcluster(P1,idx1A,P2,idx2B,eta,Nmin,iStart,jStart);
    A.M{2,1} = matrixcluster(P1,idx1B,P2,idx2A,eta,Nmin,iStart,jStart);
    A.M{2,2} = matrixcluster(P1,idx1B,P2,idx2B,eta,Nmin,iStart,jStart);
end


end