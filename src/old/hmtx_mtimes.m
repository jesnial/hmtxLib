function C = hmtx_mtimes(A,B)

% MTX_MTIMES performs the matrix multiplication. Either A or B or both can
% be H-matrices
%
% USE:
% C = hmtx_mtimes(A,B)
%
% INPUTS:
% 'A': H-matrix/fullmatrix
% 'B': H-matrix/fullmatrix
%
% OUTPUTS:
% 'C': full matrix (or H-matrix when both A and B are H-matrices) such that
%   C = A*B 
%
% NOTE:
%
% VERSION:
% Date: 04.09.2014
% Copyright(C) 2014: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:

if ishmtx(A) && isfloat(B) && ismatrix(B)
    % check input compatibility and perform product
    [nr,nc] = size(B);
    if A.ncol == nr
        C = zeros(A.nrow,nc);
        C = block_HxM(A,B,C);
    else
        error('HMTXLIB:HMTX_MTIMES','Inner matrix dimensions must agree\n');
    end
elseif isfloat(A) && ismatrix(A) && ishmtx(B)
    % check input compatibility and perform product
    [nr,nc] = size(A);
    if nc == B.nrow
        C = zeros(nr,B.ncol);
        C = block_MxH(A,B,C);
    else
        error('HMTXLIB:HMTX_MTIMES','Inner matrix dimensions must agree\n');
    end
elseif isfloat(A) && ismatrix(A) && isfloat(B) && ismatrix(B)
    C = A*B;
elseif ishmtx(A) && ishmtx(B)
    C = hmtx_mult(A,B);
end

end
