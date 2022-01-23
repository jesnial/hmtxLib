function c = hmtx_getcluster(H)

% HMTX_GETCLUSTER get the cluster IDs processing the diagonal of the
% H-matrix
%
% USE:
% c = hmtx_getcluster(H)
%
% INPUTS:
% 'H': H-matrix
%
% OUTPUTS:
% 'c': cell array with clusters
%
% NOTE:
%
% VERSION:
% Date: 22.04.2013
% Copyright(C) 2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:

% initialization
c = cell(0);

c = getcluster(H,c);

end

function c = getcluster(H,c)

% first diagonal block
if strcmpi(H.M{1,1}.type,'supermatrix')
    c = getcluster(H.M{1,1},c);
elseif any(strcmpi(H.M{1,1}.type,{'fullmatrix','rkmatrix'}))
    c{length(c)+1} = H.M{1,1}.irow;
end

% second diagonal block
if strcmpi(H.M{2,2}.type,'supermatrix')
    c = getcluster(H.M{2,2},c);
elseif any(strcmpi(H.M{2,2}.type,{'fullmatrix','rkmatrix'}))
    c{length(c)+1} = H.M{2,2}.irow;
end

end
