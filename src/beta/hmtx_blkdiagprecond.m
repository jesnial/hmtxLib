function D = hmtx_blkdiagprecond(H)

% HMTX_BLKDIAGPRECOND creates a block diagonal preconditioner by inverting
% the fullmatrices along the main diagonal
%
% USE:
% D = hmtx_blkdiagprecond(H)
%
% INPUTS:
% 'H': H-matrix structure
%
% OUTPUTS:
% 'D': supermatrix with diagonal made of the inverse of the diagonal blocks
%   of H 
%
% NOTE:
%
% VERSION:
% Date: 14.02.2014
% Copyright(C) 2014: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:

if strcmpi(H.type,'supermatrix')
    D = hmtx_create('supermatrix',H.irow,H.jcol);
    
    % create diagonal fullmatrices
    D.M{1,1} = hmtx_create('fullmatrix',H.M{1,1}.irow,H.M{1,1}.jcol);
    D.M{2,2} = hmtx_create('fullmatrix',H.M{2,2}.irow,H.M{2,2}.jcol);

    % create fake off-diagonal rkmatrices
    D.M{2,1} = hmtx_create('rkmatrix',H.M{2,1}.irow,H.M{2,1}.jcol);
    D.M{1,2} = hmtx_create('rkmatrix',H.M{1,2}.irow,H.M{1,2}.jcol);
    
    D.M{1,1} = hmtx_blkdiagprecond(H.M{1,1});
    D.M{2,2} = hmtx_blkdiagprecond(H.M{2,2});
elseif strcmpi(H.type,'fullmatrix')
    D = hmtx_create('fullmatrix',H.irow,H.jcol);
    D.M = inv(H.M);
    D.eps = H.eps;
else
    error('HMTX_BLKDIAGPRECOND: Diagonal block in rkmatrix format\n');
end

    
end

