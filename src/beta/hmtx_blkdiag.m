function D = hmtx_blkdiag(H)

% HMTX_BLKDIAG extracts the diagonal blocks
%
% USE:
% D = hmtx_blkdiag(H)
%
% INPUTS:
% 'H': H-matrix structure
%
% OUTPUTS:
% 'D': supermatrix with the diagonal blocks of H 
%
% NOTE:
%
% VERSION:
% Date: 16.02.2014
% Copyright(C) 2014: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 17.02.2014: some cleaning of the routine

if strcmpi(H.type,'supermatrix')
    D = hmtx_create('supermatrix',H.irow,H.jcol);
    
    % create fake off-diagonal rkmatrices
    D.M{2,1} = hmtx_create('rkmatrix',H.M{2,1}.irow,H.M{2,1}.jcol);
    D.M{1,2} = hmtx_create('rkmatrix',H.M{1,2}.irow,H.M{1,2}.jcol);

    % create diagonal fullmatrices
    D.M{1,1} = hmtx_blkdiag(H.M{1,1});
    D.M{2,2} = hmtx_blkdiag(H.M{2,2});
elseif strcmpi(H.type,'fullmatrix')
    D = H;
else
    error('HMTX_BLKDIAGPRECOND: Diagonal block in rkmatrix format\n');
end
    
end

