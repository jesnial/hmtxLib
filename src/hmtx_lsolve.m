function x = hmtx_lsolve(L,b)

% HMTX_LSOLVE solves a triangular system with full rhs
%
% USE:
% x = hmtx_lsolve(L,b)
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
% Date: 25.01.2013
% Copyright(C) 2013-2019: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 29.11.2019: final refurbished version
% 17.12.2019: fixed bug in the solution of the lower block

if ismatrix(b)
    % preallocation
    x = zeros(size(b));
    % recursive call
    x = lsolve(L,b,x);
else
    error('HMTX_LSOLVE: rhs must be a full vector or matrix\n');  
end

end

% effective multiplication or recursive call
function x = lsolve(L,b,x)

% Note: b and x are always the full vectors

if strcmpi(L.type,'fullmatrix')
    x(L.jcol,:) = L.M\b(L.irow,:);
elseif strcmpi(L.type,'supermatrix')
    
    % upper block solution
    x = lsolve(L.M{1,1},b,x);
    
    % rhs
    L21x = zeros(size(b));    
    L21x(L.M{2,1}.irow,:) =  hmtx_HxM(L.M{2,1},x(L.M{2,1}.jcol));
    
    % lower block solution
    x = lsolve(L.M{2,2},b-L21x,x);
end

end