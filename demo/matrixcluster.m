function A = matrixcluster(P1,idx1,P2,idx2,Nmin)

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
    A = hmtx_create('rkmatrix',idx1,idx2);
elseif length(idx1) < Nmin || length(idx2) < Nmin
    A = hmtx_create('fullmatrix',idx1,idx2);
else
    A = hmtx_create('supermatrix',idx1,idx2);
    
    % dissection
    [idx1A,idx1B] = dissect3d(P1(idx1,:));
    [idx2A,idx2B] = dissect3d(P2(idx2,:));
    
    % local-to-global
    idx1A = idx1(idx1A);
    idx1B = idx1(idx1B);
    idx2A = idx2(idx2A);
    idx2B = idx2(idx2B);
    
    % recursive call
    A.M{1,1} = matrixcluster(P1,idx1A,P2,idx2A,eta,Nmin);
    A.M{1,2} = matrixcluster(P1,idx1A,P2,idx2B,eta,Nmin);
    A.M{2,1} = matrixcluster(P1,idx1B,P2,idx2A,eta,Nmin);
    A.M{2,2} = matrixcluster(P1,idx1B,P2,idx2B,eta,Nmin);
    
end


end

