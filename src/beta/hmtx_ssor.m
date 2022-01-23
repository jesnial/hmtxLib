function [L,U] = hmtx_ssor(A,w)

% HMTX_SSOR builds the symmetric successive over-relaxation preconditioner
%
% USE:
% [L,U] = hmtx_ssor(A,w)
%
% INPUTS:
% 'A': input H-matrix
%
% OUTPUT:
% 'L': left preconditioner
% 'U': right preconditioner
%
% NOTE:
%
% VERSION:
% Date: 17.02.2014
% Copyright(C) 2014: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:

% default w
if nargin == 1
    w = 2/3;
end

% splitting A = D+L+U;
% diagonal
D = hmtx_blkdiag(A);
% lower triangular
L = hmtx_tril(A);
% upper triangular
U = hmtx_triu(A);
% inverse of diagonal
Di = hmtx_blkdiagprecond(A);

% left preconditioner  L = (D+w*L)*inv(D)
L = hmtx_muladd([],hmtx_add(D,hmtx_tH(L,w),1,1),Di);

% right preconditioner U = (D+w*L);
U = hmtx_add(D,hmtx_tH(U,w),1,1);
