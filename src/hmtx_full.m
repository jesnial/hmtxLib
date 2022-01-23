function B = hmtx_full(A)

% HMTX_FULL converts a generic H-matrix to a one-level 'fullmatrix'
%
% USE:
% B = hmtx_full(A)
%
% INPUTS:
% 'A': H-matrix
%
% OUTPUTS:
% 'B': one-level 'fullmatrix'
%
% NOTE:
%
% VERSION:
% Date: 23.12.2012
% Copyright(C) 2012-2020: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 14.01.2013: fixed bug on indices
% 21.01.2013: changed name from FULLMATRIX
% 23.01.2013: output changed into H-matrix format
% 23.01.2013: uses HMTX_CREATE to create H-matrix structure
% 24.01.2013: changed call to ASSEMBLEFULL
% 24.01.2013: added 'eps' calculation
% 05.01.2020: uses HMTX_CHANGEFORMAT

% create output structure
B = hmtx_create('fullmatrix',A.irow,A.jcol);

% change format
B = hmtx_changeformat(A,B);

% B.M = zeros(A.nrow,A.ncol);
% B = assemblefull(A,B);

end

% effective loading of full matrix or recursive call
function F = assemblefull(A,F)

if strcmpi(A.type,'supermatrix')
    F = assemblefull(A.M{1,1},F);
    F = assemblefull(A.M{1,2},F);
    F = assemblefull(A.M{2,1},F);
    F = assemblefull(A.M{2,2},F);
elseif strcmpi(A.type,'fullmatrix')
    % get local numbering
    [~,irow] = ismember(A.irow,F.irow);
    [~,jcol] = ismember(A.jcol,F.jcol);
    % load values
    F.M(irow,jcol) = A.M;
elseif strcmpi(A.type,'rkmatrix')
    % get local numbering
    [~,irow] = ismember(A.irow,F.irow);
    [~,jcol] = ismember(A.jcol,F.jcol);
    % load values
    F.M(irow,jcol) = A.U*A.V';
end
F.eps = max([F.eps,A.eps]);

end
