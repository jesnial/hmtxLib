function C = block_HxM(H,M,C)

% BLOCK_HxM performs the H-matrix times fullmatrix product
%
% USE:
% C = block_HxM(H,M,C)
%
% INPUTS:
% 'H': H-matrix structure
% 'M': full matrix
% 'C': output fullmatrix
%
% OUTPUTS:
% 'C': full matrix A = H*M
%
% NOTE:
%
% VERSION:
% Date: 21.12.2012
% Copyright(C) 2012-2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 08.01.2013: grouped partition of vector
% 10.01.2013: fixed bug of indexing
% 17.01.2013: generalized for matrices
% 21.01.2013: changed name from HMTX_AX
% 21.01.2013: changed oreder of operations in the case of 'rkmatrix'

if strcmpi(H.type,'supermatrix')
    C = block_HxM(H.M{1,1},M,C);
    C = block_HxM(H.M{1,2},M,C);
    C = block_HxM(H.M{2,1},M,C);
    C = block_HxM(H.M{2,2},M,C);
elseif strcmpi(H.type,'fullmatrix') && ~isempty(H.M)
    C(H.irow,:) = C(H.irow,:)+H.M*M(H.jcol,:);
elseif strcmpi(H.type,'rkmatrix') && ~isempty(H.U)
    C(H.irow,:) = C(H.irow,:)+H.U*(M(H.jcol,:)'*H.V)';
end

end