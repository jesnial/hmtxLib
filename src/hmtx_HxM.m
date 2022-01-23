function A = hmtx_HxM(H,M)

% HMTX_HxM performs the H-matrix times double matrix product
%
% USE:
% A = hmtx_HxM(H,M)
%
% INPUTS:
% 'H': H-matrix structure
% 'M': full matrix
%
% OUTPUTS:
% 'A': full matrix A = H*M
%
% NOTE:
%
% VERSION:
% Date: 21.12.2012
% Copyright(C) 2012-2019: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 08.01.2013: grouped partition of vector
% 10.01.2013: fixed bug of indexing
% 17.01.2013: generalized for matrices
% 21.01.2013: changed name from HMTX_AX
% 21.01.2013: changed oreder of operations in the case of 'rkmatrix'
% 29.11.2019: generalized for supermatrices with partial indices

A = zeros(H.nrow,size(M,2));
irowGlo = H.irow;
jcolGlo = H.jcol;
A = block_HxM(H,M,A,irowGlo,jcolGlo);

end

% effective multiplication or recursive call
function A = block_HxM(H,M,A,irowGlo,jcolGlo)

if strcmpi(H.type,'supermatrix')
    A = block_HxM(H.M{1,1},M,A,irowGlo,jcolGlo);
    A = block_HxM(H.M{1,2},M,A,irowGlo,jcolGlo);
    A = block_HxM(H.M{2,1},M,A,irowGlo,jcolGlo);
    A = block_HxM(H.M{2,2},M,A,irowGlo,jcolGlo);
elseif strcmpi(H.type,'fullmatrix') && ~isempty(H.M)
    [~,irowLoc] = ismember(H.irow,irowGlo);
    [~,jcolLoc] = ismember(H.jcol,jcolGlo);
    A(irowLoc,:) = A(irowLoc,:)+H.M*M(jcolLoc,:);
elseif strcmpi(H.type,'rkmatrix') && ~isempty(H.U)
    [~,irowLoc] = ismember(H.irow,irowGlo);
    [~,jcolLoc] = ismember(H.jcol,jcolGlo);
    A(irowLoc,:) = A(irowLoc,:)+H.U*(M(jcolLoc,:)'*H.V)';
end

end