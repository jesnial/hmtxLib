function [H,b] = hmtx_fixedboundary(H,b,idfix,xfix)

% HMTX_FIXEDBOUNDARY nullifies the rows specified and put 1s on the diagonal
%
% USE:
% [H,b] = hmtx_fixedboundary(H,b,idfix,ixfix)
%
% INPUTS:
% 'H': H-matrix
% 'b': rhs
% 'idfix': indices of elems with assigned values (e.g. Dirichlet boundary conditions)
% 'xfix': values to be assigned to elems
% 
% OUTPUTS:
% 'H': modified H-matrix
% 'b': modified rhs
%
% NOTE:
%
% VERSION:
% Date: 22.04.2013
% Copyright(C) 2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
%

% rhs
b(idfix) = xfix;

% matrix
H = fixedboundary(H,idfix);

end

% recursive function
function H = fixedboundary(H,idfix)

if strcmpi(H.type,'supermatrix')
    % recursive call
    H.M{1,1} = fixedboundary(H.M{1,1},idfix);
    H.M{1,2} = fixedboundary(H.M{1,2},idfix);
    H.M{2,1} = fixedboundary(H.M{2,1},idfix);
    H.M{2,2} = fixedboundary(H.M{2,2},idfix);
elseif strcmpi(H.type,'fullmatrix')
    % find rows in local numbering (unsorted)
    ifix = find(ismember(H.irow,idfix));
    % nullify row
    H.M(ifix,:) = 0;
    % put 1s on the diagonal
    [tf, loc] = ismember(H.irow(ifix),H.jcol);
    % find diagonal elements in local numbering (no ordering assumed)
    idiag = ifix(tf);
    jdiag = loc(loc ~= 0);
    H.M(sub2ind([H.nrow,H.ncol],idiag,jdiag)) = 1;
elseif strcmpi(H.type,'rkmatrix')
    % find rows in local numbering (unsorted)
    ifix = ismember(H.irow,idfix);
    % nullify row
    H.U(ifix,:) = 0;
end

end