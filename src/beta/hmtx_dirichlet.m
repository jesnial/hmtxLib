function [H,b] = hmtx_dirichlet(H,b,idfix,xfix)

% HMTX_DIRICHLET modifies the coefficient matrix end rhs imposing the value
%   of some unknowns
%
% USE:
% [H,b] = hmtx_dirichlet(H,b,idfix,ixfix)
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
% Algorithm: assuming to set the value of x(i) to k
% 1) the i-th row is nullified
% 2) the diagonal value is set to 1
% 3) the i-th entry of the rhs is set to the known value
% 4) the i-th column are nullified
% 5) the rhs (with the exception of the i-th entry) is decreased of H(:,i)*k
%
% VERSION:
% Date: 05.03.2014
% Copyright(C) 2014: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:

% matrix
[H,b] = dirichlet(H,b,int32(idfix(:)),xfix(:));

% rhs
b(idfix) = xfix;

end

% recursive function
function [H,b] = dirichlet(H,b,idfix,xfix)

if strcmpi(H.type,'supermatrix')
    % recursive call
    [H.M{1,1},b] = dirichlet(H.M{1,1},b,idfix,xfix);
    [H.M{1,2},b] = dirichlet(H.M{1,2},b,idfix,xfix);
    [H.M{2,1},b] = dirichlet(H.M{2,1},b,idfix,xfix);
    [H.M{2,2},b] = dirichlet(H.M{2,2},b,idfix,xfix);
elseif strcmpi(H.type,'fullmatrix')
    % find rows in local numbering (unsorted)
    ifix = find(ismembc(H.irow,idfix));
    % find cols in local numbering (unsorted)
    loc = ismembc2(H.jcol,idfix);
    jfix = find(loc > 0);
    
    % modify rhs
    if ~isempty(jfix)
        b(H.irow) = b(H.irow)-H.M(:,jfix)*xfix(loc(loc > 0));
    end
    
    % nullify row
    H.M(ifix,:) = 0;
    % nullify col
    H.M(:,jfix) = 0;
    
    % put 1s on the diagonal
    loc = ismembc2(H.irow(ifix),H.jcol);
    % find diagonal elements in local numbering (no ordering assumed)
    jdiag = loc(loc ~= 0);
    idiag = ifix(loc ~= 0);
    H.M(sub2ind([H.nrow,H.ncol],idiag,jdiag)) = 1;
elseif strcmpi(H.type,'rkmatrix')
    % find rows in local numbering (unsorted)
    ifix = ismembc(H.irow,idfix);
    % find cols in local numbering (unsorted)
    loc = ismembc2(H.jcol,idfix);
    jfix = find(loc > 0);
    
    % modify rhs
    if ~isempty(jfix)
        b(H.irow) = b(H.irow)-H.U*H.V(jfix,:)'*xfix(loc(loc > 0));
    end
    
    % nullify row
    H.U(ifix,:) = 0;
    % nullify col
    H.V(jfix,:) = 0;
end

end