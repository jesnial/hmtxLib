function At = hmtx_ctranspose(A)

% HMTX_CTRANSPOSE calculates the complex conjugate transpose of an H-matrix
%
% USE:
% At = hmtx_ctranspose(A)
%
% INPUTS:
% 'A': H-matrix structure, as created by HMTX_CLUSTER and filled by HMTX_FILL
%
% OUTPUTS:
% 'At': complex conjugate transpose H-matrix of A
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
    At = hmtx_create('supermatrix',A.jcol,A.irow);
    At.eps = A.eps;
    
    % get transpose of diagonal blocks
    At.M{1,1} = hmtx_ctranspose(A.M{1,1});
    At.M{2,2} = hmtx_ctranspose(A.M{2,2});
    
    % exchange other blocks
    At.M{1,2} = hmtx_ctranspose(A.M{2,1});
    At.M{2,1} = hmtx_ctranspose(A.M{1,2});
    
elseif strcmpi(A.type,'fullmatrix')
    At = hmtx_create('fullmatrix',A.jcol,A.irow);
    At.M = A.M';
    At.eps = A.eps;
    
elseif strcmpi(A.type,'rkmatrix')
    At = hmtx_create('rkmatrix',A.jcol,A.irow);
    At.k = A.k;
    At.kMax = A.kMax;
    At.eps = A.eps;
    At.U = A.V;
    At.V = A.U;
end

end