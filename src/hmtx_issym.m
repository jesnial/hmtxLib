function Bool = hmtx_issym(A, tol)

% HMTX_ISSYM checks if an H-matrix is symmetric 
% 
% USE:
% Bool = hmtx_issym(A, tol)
%
% INPUTS:
% 'A': H-matrix structure, as created by HMTX_CLUSTER and filled by HMTX_FILL
% 'tol' : tolerance of the method
%
% OUTPUTS:
% 'Bool': true if A is symmetric, false otherwise
%
% NOTE:
% For the Cholesky decomposition, one should also check if A is positive
% definite
%
% VERSION:
% Date: 15.01.2020
% Copyright(C) 2020: Jessie Levillain (jessielevillain@gmail.com)
%
% HISTORY:

AT = hmtx_transpose(A);
AAT = getfield(hmtx_full(hmtx_add(A, AT, '-')),'M');
relnorm = norm(AAT,'fro');
if relnorm < tol
    Bool = true ;
else
    Bool = false;
end
end
