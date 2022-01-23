function B = hmtx_makesparse(A)
    
% HMTX_MAKESPARSE converts a generic H-matrix to a sparse matrix using only
% the fullblocks
%
% USE:
% B = hmtx_makesparse(A)
%
% INPUTS:
% 'A': H-matrix
%
% OUTPUTS:
% 'B': sparse matrix made of full blocks
%
% NOTE:
%
% VERSION:
% Date: 27.12.2019
% Copyright(C) 2019: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:

irow = [];
jcol = [];
b = [];

% recursive call
[irow,jcol,b] = getfullblocks(A,irow,jcol,b);

% create sparse matrix
B = sparse(double(irow),double(jcol),b,A.nrow,A.ncol);

end

% effective loading of full matrix or recursive call
function [irow,jcol,b] = getfullblocks(A,irow,jcol,b)

if strcmpi(A.type,'supermatrix')
    [irow,jcol,b] = getfullblocks(A.M{1,1},irow,jcol,b);
    [irow,jcol,b] = getfullblocks(A.M{1,2},irow,jcol,b);
    [irow,jcol,b] = getfullblocks(A.M{2,1},irow,jcol,b);
    [irow,jcol,b] = getfullblocks(A.M{2,2},irow,jcol,b);
elseif strcmpi(A.type,'fullmatrix')
    % load indices and values
    [irowloc,jcolloc] = meshgrid(A.irow,A.jcol);
    irow = [irow; irowloc(:)];
    jcol = [jcol; jcolloc(:)];
    b = [b; A.M(:)];
end

end
