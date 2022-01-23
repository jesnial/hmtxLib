function B = hmtx_copystruct(A)

% HMTX_COPYSTRUCT creates an empty matrix with the same structure of A
%
% USE:
% B = hmtx_copystruct(A)
%
% INPUTS:
% 'A': H-matrix structure, as created by HMTX_CLUSTER
%
% OUTPUTS:
% 'B': H-matrix with the same structure of A
%
% NOTE:
% Recursive function
%
% VERSION:
% Date: 06.02.2014
% Copyright(C) 2014-2020: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 04.09.2014: help added
% 03.01.2020: added KMAX field for rkmatrix;

if strcmpi(A.type,'supermatrix')
    B = hmtx_create('supermatrix',A.irow,A.jcol);
    B.eps = A.eps;
    % copy blocks
    B.M{1,1} = hmtx_copystruct(A.M{1,1});
    B.M{2,1} = hmtx_copystruct(A.M{2,1});
    B.M{1,2} = hmtx_copystruct(A.M{1,2});
    B.M{2,2} = hmtx_copystruct(A.M{2,2});
elseif strcmpi(A.type,'fullmatrix')
    B = hmtx_create('fullmatrix',A.irow,A.jcol);
    B.eps = A.eps;
elseif strcmpi(A.type,'rkmatrix')
    B = hmtx_create('rkmatrix',A.irow,A.jcol);
    B.eps = A.eps;
    B.kMax = A.kMax;
end