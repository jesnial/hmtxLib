function C = block_MxH(M,H,C)

% BLOCK_MxH performs the fullmatrix times H-matrix product
%
% USE:
% C = hmtx_MxH(M,H,C)
%
% INPUTS:
% 'M': full matrix
% 'H': H-matrix structure
% 'C': output fullmatrix
%
% OUTPUTS:
% 'C': full matrix A = M*H
%
% NOTE:
%
% VERSION:
% Date: 21.01.2013
% Copyright(C) 2012-2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:

if strcmpi(H.type,'supermatrix')
    C = block_MxH(M,H.M{1,1},C);
    C = block_MxH(M,H.M{2,1},C);
    C = block_MxH(M,H.M{1,2},C);
    C = block_MxH(M,H.M{2,2},C);
elseif strcmpi(H.type,'fullmatrix')
    C(:,H.jcol) = C(:,H.jcol)+M(:,H.irow)*H.M;
elseif strcmpi(H.type,'rkmatrix')
    C(:,H.jcol) = C(:,H.jcol)+(M(:,H.irow)*H.U)*H.V';
end

end