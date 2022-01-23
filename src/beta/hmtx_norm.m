function n = hmtx_norm(H)

% HMTX_NORM calculates the Frobenius norm of an H-matrix
%
% USE:
% n = hmtx_norm(H)
%
% INPUTS:
% 'A': H-matrix
%
% OUTPUTS:
% 'n': Frobenius norm
%
% VERSION:
% Date: 17.02.2014
% Copyright(C) 2014: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:

% evaluate square of frobenius norm
n2 = fronorm(H);
n = sqrt(n2);

end

function n2 = fronorm(H)

if strcmpi(H.type,'supermatrix')
    n11 = fronorm(H.M{1,1});
    n21 = fronorm(H.M{2,1});
    n12 = fronorm(H.M{1,2});
    n22 = fronorm(H.M{2,2});
    n2 = n11+n21+n12+n22;
elseif strcmpi(H.type,'fullmatrix')
    n2 = sum(H.M(:).^2);
elseif strcmpi(H.type,'rkmatrix')
    if H.k ~= 0
        M = H.U*H.V';
        n2 = sum(M(:).^2);
    else
        n2 = 0;
    end
end


end