function x = hmtx_ltsolve(L,b)

% HMTX_LTSOLVE solves a triangular system with full rhs such as x = L'\b
%
% USE:
% x = hmtx_ltsolve(L,b)
%
% INPUTS:
% 'L': lower triangular H-matrix
% 'b': full rhs ('b' can be a matrix)
%
% OUTPUTS:
% 'x': solution of the triangular system
%
% NOTE:
%
% VERSION:
% Date: 17.01.2020
% Copyright(C) 2020: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:

if ismatrix(b)
    % preallocation
    x = zeros(size(b));
    % recursive call
    x = ltsolve(L,b,x);
else
    error('HMTX_LSOLVE: rhs must be a full vector or matrix\n');  
end

end

% effective multiplication or recursive call
function x = ltsolve(L,b,x)

% Note: b and x are always the full vectors

if strcmpi(L.type,'fullmatrix')
    x(L.irow,:) = L.M'\b(L.jcol,:);
elseif strcmpi(L.type,'supermatrix')
    
    % upper block solution
    x = ltsolve(L.M{2,2},b,x);
    
    % rhs
    Lt21x = zeros(size(b));    
    Lt21x(L.M{2,1}.jcol,:) =  hmtx_HxM(hmtx_herm(L.M{2,1}),x(L.M{2,1}.irow));
    
    % lower block solution
    x = ltsolve(L.M{1,1},b-Lt21x,x);
end

end