function A = hmtx_MxH(M,H)

% HMTX_MxH performs the fullmatrix times H-matrix product 
%
% USE:
% A = hmtx_MxH(M,H)
%
% INPUTS:
% 'M': full matrix
% 'H': H-matrix structure
%
% OUTPUTS:
% 'A': full matrix A = M*H
%
% NOTE:
%
% VERSION:
% Date: 21.01.2013
% Copyright(C) 2012-2019: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 29.11.2019: generalized for supermatrices with partial indices

A = zeros(size(M,1),H.ncol);
irowGlo = H.irow;
jcolGlo = H.jcol;
A = block_MxH(M,H,A,irowGlo,jcolGlo);

end

% effective multiplication or recursive call
function A = block_MxH(M,H,A,irowGlo,jcolGlo)

if strcmpi(H.type,'supermatrix')
    A = block_MxH(M,H.M{1,1},A,irowGlo,jcolGlo);
    A = block_MxH(M,H.M{2,1},A,irowGlo,jcolGlo);
    A = block_MxH(M,H.M{1,2},A,irowGlo,jcolGlo);
    A = block_MxH(M,H.M{2,2},A,irowGlo,jcolGlo);
elseif strcmpi(H.type,'fullmatrix') && ~isempty(H.M)
    [~,irowLoc] = ismember(H.irow,irowGlo);
    [~,jcolLoc] = ismember(H.jcol,jcolGlo);
    A(:,jcolLoc) = A(:,jcolLoc)+M(:,irowLoc)*H.M;
elseif strcmpi(H.type,'rkmatrix') && ~isempty(H.U)
    [~,irowLoc] = ismember(H.irow,irowGlo);
    [~,jcolLoc] = ismember(H.jcol,jcolGlo);
    A(:,jcolLoc) = A(:,jcolLoc)+(M(:,irowLoc)*H.U)*H.V';
end

end