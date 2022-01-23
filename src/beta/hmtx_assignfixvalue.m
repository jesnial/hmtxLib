function [H,b,idunk] = hmtx_assignfixvalue(H,b,idfix,xfix)

% HMTX_ASSIGNFIXVALUE assigns to certain idfix a prefixed xfix and rearrange
%                the solving matrix and rhs taking into account equivalues
%
% USE:
% [H,b,idunk] = hmtx_assignfixvalue(H,b,idfix,xfix)
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
% 'idunk': absolute indices of unknowns
%
% NOTE:
% 
%
% VERSION:
% Date: 06.11.2015
% Copyright(C) 2015: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
%

% remove repeated assignments
[idfix,ielems] = unique(idfix);

% fixed values
if ~isempty(idfix)
    % create complete vector with only xfix as nonzeros
    x = zeros(size(b));
    x(idfix) = xfix;
    %b = b-repmat(A(:,idfix)*xfix(ielems),1,size(b,2));
    b = b-hmtx_HxM(H,x);
end

% rhs
b(idfix) = xfix;

% matrix
H = deleteRowsCols(H,idfix);

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