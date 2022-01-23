function C = hmtx_add(A,B,op)

% HMTX_ADD sums two H-matrices with the same structure
%
% USE:
% C = hmtx_add(A,B)
%
% INPUTS:
% 'A': H-matrix, as created by HMTX_CLUSTER and filled by HMTX_FILL
% 'B': H-matrix, as created by HMTX_CLUSTER and filled by HMTX_FILL
% 'op': operator: '+'/'-'
%
% OUTPUTS:
% 'C': H-matrix such that C = A+B or C = A-B according to 'op'
%
% NOTE:
% The result of the sum of two rkmatrices is truncated by using
% rSVD_rkmatrix to keep the max EPS and/or min KMAX
%
% VERSION:
% Date: 14.01.2013
% Copyright(C) 2013-2020: Fabio Freschi (fabio.freschi@polito.it)
%                         Jessie Levillain (jessielevillain@gmail.com)
%
% HISTORY:
% 22.01.2013: added case of type mismatch
% 22.01.2013: added OP as input to take into account the subtraction
% 22.01.2013: correction of EPS assignment
% 23.01.2013: uses HMTX_CREATE to create H-matrix structure
% 13.11.2019: sign correction for Rk matrices
% 16.12.2019: new use of FULL2SUPER and RKMATRIX2SUPER
% 02.01.2020: added k and eps fields to rkmatrix
% 03.01.2020: added kMax-truncated addition
% 05.01.2020: reoved unused cases

if strcmpi(op,'+')
    t = 1;
elseif strcmpi(op,'-')
    t = -1;
end

if strcmpi(A.type,'supermatrix') && strcmpi(B.type,'supermatrix')
    C = hmtx_create('supermatrix',A.irow,A.jcol);
    C.M{1,1} = hmtx_add(A.M{1,1},B.M{1,1},op);
    C.M{1,2} = hmtx_add(A.M{1,2},B.M{1,2},op);
    C.M{2,1} = hmtx_add(A.M{2,1},B.M{2,1},op);
    C.M{2,2} = hmtx_add(A.M{2,2},B.M{2,2},op);
    C.eps = max([A.eps B.eps]);
elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'fullmatrix')
    C = hmtx_create('fullmatrix',A.irow,A.jcol);
    C.M = A.M+t*B.M;
    C.eps = max([A.eps B.eps]);
elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'rkmatrix')
    C = hmtx_create('rkmatrix',A.irow,A.jcol);
    C.U = [A.U t*B.U];
    C.V = [A.V B.V];
    C.k = A.k+B.k;
    eps = max([A.eps B.eps]);
    kMax = min([A.kMax B.kMax]);
    C = rSVD_rkmatrix(C,eps,kMax);
else
    error('Matlab:hmtx_add is called with %s %s matrix blocks\n',A.type,B.type);
end
% elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'fullmatrix')
%     C = hmtx_create('fullmatrix',A.irow,A.jcol);
%     C.M = A.U*A.V'+t*B.M;
%     C.eps = A.eps;
% elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'rkmatrix')
%     C = hmtx_create('fullmatrix',A.irow,A.jcol);
%     C.M = A.M+t*B.U*B.V';
%     C.eps = B.eps;
% elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'rkmatrix')
%     B2 = hmtx_copystruct(A);
%     B2 = rkmatrix2super(B,B2);
%     C = hmtx_add(A,B2,op);
% elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'supermatrix')
%     A2 = hmtx_copystruct(B);
%     A2 = rkmatrix2super(A,A2);
%     C = hmtx_add(A2,B,op);
% elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'fullmatrix')
%     B2 = hmtx_copystruct(A);
%     B2 = full2super(B,B2);
%     C = hmtx_add(A,B2,op);
% elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'supermatrix')
%     A2 = hmtx_copystruct(B);
%     A2 = full2super(A,A2);
%     C = hmtx_add(A2,B,op);
% end

end