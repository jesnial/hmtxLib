function L = hmtx_chol(A)

% Calculates an (approximate) cholesky decomposition of an H-matrix A such
% that A ~ L*L'
%
% USE:
% L = hmtx_chol(A)
%
% INPUTS:
% 'A': symmetric H-matrix
%
% OUTPUTS:
% 'L': lower triangular cholesky factor
%
% NOTE:
% This function does not make any check on the matrix and it works on the
% lower triangular part
% See M. Bebendorf, 'Hierarchical Matrices', p. 87
%
% VERSION:
% Date: 15.01.2020
% Copyright(C) 2020: Jessie Levillain (jessielevillain@gmail.com)
%                    Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 17.01.2020: removed 'flag' from input list

if strcmpi(A.type,'supermatrix')
     L = hmtx_create('supermatrix', A.irow, A.jcol);
     L.eps = A.eps;
     
     % create blocks full of 0 to avoid problems
     L.M{1,2} = hmtx_create('rkmatrix', A.M{1,2}.irow, A.M{1,2}.jcol);
     L.M{1,2}.eps = 0;
     L.M{1,2}.k = 1;
     L.M{1,2}.U = zeros(A.M{1,2}.nrow, L.M{1,2}.k);
     L.M{1,2}.V = zeros(A.M{1,2}.ncol, L.M{1,2}.k);
     
     % find LU for first block
     [L.M{1,1}] = hmtx_chol(A.M{1,1});
     
     % solve for diagonal block
     L.M{2,1} = hmtx_ctranspose(hmtx_leftsolve(L.M{1,1}, A.M{1,2}));
     
     % find LL' for last blocks (L22U22 = A22 - L21*L12' = B)
     B = hmtx_copystruct(A.M{2,2});
     B = hmtx_mult(L.M{2,1},hmtx_ctranspose(L.M{2,1}),B);
     L.M{2,2} = hmtx_chol(hmtx_add(A.M{2,2},B,'-')); 
elseif strcmpi(A.type,'fullmatrix')
    % use standard chol decomposition
    L = hmtx_create('fullmatrix',A.irow,A.jcol);
    [L.M,flag] = chol(A.M, 'lower');
    if flag > 0
        error('hmtx_cholesky:A diagonal H-matrix block is not symmetric positive definite\n')
    end
end

    


end

