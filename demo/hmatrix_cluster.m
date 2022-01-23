function A = hmatrix_cluster(P1,varargin)

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


% check varargin
for i = 1:2:length(varargin)
    vin = varargin{i};
    if strcmpi(vin,'eta')
        eta = varargin{i+1};
    elseif strcmpi(vin,'Nmin')
        Nmin = varargin{i+1};
    elseif strcmpi(vin,'P2')
        P2 = varargin{i+1};
    else
        warning('MATLAB:hmtx_cluster','Wrong VARARGIN parameter, skip ''%s = %s''\n',vin,num2str(varargin{i+1}));
    end
end

% create matrix cluster structure
A = matrixcluster(P1,(1:size(P1,1)).',P2,(1:size(P2,1)).',eta,Nmin);

end

function C = matrixcluster(P1,idx1,P2,idx2,eta,Nmin)

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

% get cluster sizes
[N,~] = size(idx1);
[M,~] = size(idx2);

% check admissibility of blocks
if clusteradmissibility(P1(idx1,:),P2(idx2,:),eta)
    C = hmatrix;
    C.sz = [N M];
    C.admissible = true;
elseif length(idx1) < Nmin || length(idx2) < Nmin
    C = hmatrix;
    C.sz = [N M];
    C.admissible = false;
else
    C = hmatrix;
    C.sz = [N M];
    C.admissible = false;
    
    % dissection
    [idx1A,idx1B] = cluster_split(P1(idx1,:)); 
    [idx2A,idx2B] = cluster_split(P2(idx2,:));
    
%     % local-to-global
%     idx1A = idx1(idx1A);
%     idx1B = idx1(idx1B); 
%     idx2A = idx2(idx2A);
%     idx2B = idx2(idx2B);
    
    % recursive call
    C.A11 = matrixcluster(P1,idx1A.',P2,idx2A.',eta,Nmin);
    C.A12 = matrixcluster(P1,idx1A.',P2,idx2B.',eta,Nmin);
    C.A21 = matrixcluster(P1,idx1B.',P2,idx2A.',eta,Nmin);
    C.A22 = matrixcluster(P1,idx1B.',P2,idx2B.',eta,Nmin);
    
end


end

function [idxA,idxB] = cluster_split(P)

[Nidx,~] = size(P); %get cluster size
    C_center=sum(P,2)/Nidx;
    C_radii =sqrt((P(:,1)-C_center(1)).^2+...
                  (P(:,2)-C_center(2)).^2+...
                  (P(:,3)-C_center(3)).^2);
 %cluster radius

N1 = sum(P<C_radii+C_center) %get number of points at which to split

% splitting indices 
idxA = 1:N1;
idxB = N1+1:Nidx;

end
