function x = hmtx_usolve(U,b)

% HMTX_USOLVE solves a triangular system with full rhs
%
% USE:
% x = hmtx_usolve(U,b)
%
% INPUTS:
% 'U': upper triangular H-matrix
% 'b': full rhs ('b' can be a matrix)
%
% OUTPUTS:
% 'x': solution of the triangular system
%
% NOTE:
%
% VERSION:
% Date: 17.02.2014
% Copyright(C) 2014-2019: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 24.02.2014: uses additional coefficients in HMTX_MVM
% 29.11.2019: final refurbished version
% 17.12.2019: fixed bug in the solution of the lower block

if ismatrix(b)
    % preallocation
    x = zeros(size(b));
    % recursive call
    x = usolve(U,b,x);
else
    error('HMTX_USOLVE: rhs must be a full vector or matrix\n');  
end

end

function x = usolve(U,b,x)

% Note: b and x are always the full vectors

if strcmpi(U.type,'fullmatrix')
    x(U.jcol,:) = U.M\b(U.irow,:);
elseif strcmpi(U.type,'supermatrix')
    
    % lower block solution
    x = usolve(U.M{2,2},b,x);
    
    % rhs
    U12x = zeros(size(b));
    U12x(U.M{1,2}.irow,:) = hmtx_HxM(U.M{1,2},x(U.M{1,2}.jcol));
    
    % upper block solution
    x = usolve(U.M{1,1},b-U12x,x);
end

end

