function [qr,qc] = hmtx_getclustertree(H)

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
qr = zeros(H.nrow,1);
% column index
qc = zeros(H.ncol,1);

% counters
ir = 1;
ic = 1;

% call recursive function
[qr,qc] = getclustertree(H,qr,qc,ir,ic);

end

function [qr,qc,ir,ic] = getclustertree(H,qr,qc,ir,ic)

if strcmpi(H.type,'supermatrix')
    [qr,qc,ir,ic] = getclustertree(H.M{1,1},qr,qc,ir,ic);
    [qr,qc,ir,ic] = getclustertree(H.M{2,2},qr,qc,ir,ic);
elseif strcmpi(H.type,'fullmatrix') || strcmpi(H.type,'rkmatrix')
    % update row
    qr(H.irow) = ir;
    ir = ir+1;
    % update col
    qc(H.jcol) = ic;
    ic = ic+1;
end

end