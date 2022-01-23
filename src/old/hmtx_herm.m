function AH = hmtx_ctranspose(A)

% HMTX_CTRANSPOSE calculates the 
%
% USE:
% AH = hmtx_ctranspose(A)
%
% INPUTS:
% 'A': H-matrix structure, as created by HMTX_CLUSTER and filled by HMTX_FILL
%
% OUTPUTS:
% 'AH': hermitian H-matrix of A
%
% NOTE:
% Recursive function
%
% VERSION:
% Date: 15.01.2020
% Copyright(C) 2020: Jessie Levillain (jessielevillain@gmail.com)
%                    Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 17.01.2020: added KMAX 

if strcmpi(A.type,'supermatrix')
    AH = hmtx_create('supermatrix',A.jcol,A.irow);
    AH.eps = A.eps;
    
    % get transpose of diagonal blocks
    AH.M{1,1} = hmtx_ctranspose(A.M{1,1});
    AH.M{2,2} = hmtx_ctranspose(A.M{2,2});
    
    % exchange other blocks
    AH.M{1,2} = hmtx_ctranspose(A.M{2,1});
    AH.M{2,1} = hmtx_ctranspose(A.M{1,2});
    
elseif strcmpi(A.type,'fullmatrix')
    AH = hmtx_create('fullmatrix',A.jcol,A.irow);
    AH.M = A.M';
    AH.eps = A.eps;
    
elseif strcmpi(A.type,'rkmatrix')
    AH = hmtx_create('rkmatrix',A.jcol,A.irow);
    AH.k = A.k;
    AH.kMax = A.kMax;
    AH.eps = A.eps;
    AH.U = A.V;
    AH.V = A.U;
end

end