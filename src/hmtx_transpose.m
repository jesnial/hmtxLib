function AT = hmtx_transpose(A)

% HMTX_TRANSPOSE transposes (non conjugate transpose) an H-matrix
% 
% USE:
% AT = hmtx_transpose(A)
%
% INPUTS:
% 'A': H-matrix structure, as created by HMTX_CLUSTER and filled by HMTX_FILL
%
% OUTPUTS:
% 'AT': transposed H-matrix of A
%
% NOTE:
% Recursive function
%
% VERSION:
% Date: 29.11.2019
% Copyright(C) 2019-2020 : Jessie Levillain (jessielevillain@gmail.com)
%                          Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 13.12.2019: fixed eps of fullmatrix
% 15.01.2020: fixed transposition for complex rk blocks
% 17.01.2020: added KMAX 

if strcmpi(A.type,'supermatrix')
     AT = hmtx_create('supermatrix',A.jcol,A.irow);
     AT.eps = A.eps;

     % get transpose of diagonal blocks
     AT.M{1,1} = hmtx_transpose(A.M{1,1});
     AT.M{2,2} = hmtx_transpose(A.M{2,2});
     
     % exchange other blocks 
     AT.M{1,2} = hmtx_transpose(A.M{2,1});
     AT.M{2,1} = hmtx_transpose(A.M{1,2});
     
elseif strcmpi(A.type,'fullmatrix')
     AT = hmtx_create('fullmatrix',A.jcol,A.irow);
     AT.M = transpose(A.M);
     AT.eps = A.eps;
     
elseif strcmpi(A.type,'rkmatrix')
    AT = hmtx_create('rkmatrix',A.jcol,A.irow);
    AT.k = A.k;
    AT.kMax = A.kMax;
    AT.eps = A.eps;
    AT.U = transpose(A.V');
    AT.V = transpose(A.U');
end

end