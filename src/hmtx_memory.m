function [nB,nE] = hmtx_memory(H)

% HMTX_MEMORY returns the bytes and number of entries of an H-matrix
%
% USE:
% [nB,nE] = hmtx_memory(H)
%
% INPUTS:
% 'H': H-matrix
%
% OUTPUTS:
% 'nB': gross memory in bytes (including overhead)
% 'nE': number of entries
%
% NOTE:
%
% VERSION:
% Date: 16.04.2013
% Copyright(C) 2013-2019: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 02.12.2019: name changed from HMTX_MBYTE
% 02.12.2019: added number of entries as output

% get gross memory (including overhead)
nB = getfield(whos('H'),'bytes');

% initialization
nE = 0;

% recursive function
nE = numberOfEntries(H,nE);

end

function nE = numberOfEntries(H,nE)

if strcmpi(H.type,'supermatrix')
    nE = numberOfEntries(H.M{1,1},nE);
    nE = numberOfEntries(H.M{1,2},nE);
    nE = numberOfEntries(H.M{2,1},nE);
    nE = numberOfEntries(H.M{2,2},nE);
elseif strcmpi(H.type,'fullmatrix')
    nE = nE+H.nrow*H.ncol;
elseif strcmpi(H.type,'rkmatrix')
    nE = nE+(H.nrow+H.ncol)*H.k;
end

end
