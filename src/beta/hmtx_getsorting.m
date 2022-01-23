function [sr,sc] = hmtx_getsorting(H)

% HMTX_GETCLUSTERTREE gets the cluster tree along rows and columns of an H
% matrix
%
% USE:
% [qr,qc] = hmtx_getclustertree(H)
%
% INPUTS:
% 'H': H matrix
%
% OUTPUTS:
% 'qr': vector of size H.nrow with scalar indices of the row cluster
% 'qc': vector of size H.ncol with scalar indices of the col cluster
%
% NOTE:
%
% VERSION:
% Date: 02.12.2019
% Copyright(C) 2019: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:

% row index
sr = [];
% column index
sc = [];

% call recursive function
[sr,sc] = getclustertree(H,sr,sc);

end

function [sr,sc] = getclustertree(H,sr,sc)

if strcmpi(H.type,'supermatrix')
    [sr,sc] = getclustertree(H.M{1,1},sr,sc);
    [sr,sc] = getclustertree(H.M{2,2},sr,sc);
elseif strcmpi(H.type,'fullmatrix') || strcmpi(H.type,'rkmatrix')
    % update row
    sr = [sr; H.irow];
    % update col
    sc = [sc; H.jcol];
end

end