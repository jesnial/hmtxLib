function C = hmtx_mult(A,B)

% C = A*B

% 21.01.2013

% 23.01.2013: uses HMTX_CREATE to create H-matrix structure

% check dimenisons
%
% HISTORY :
%
% 29.11.2019 : fixed cross-type cases. supermatrix and rk becomes rk
% 03.12.2019 : added recompression


if A.ncol ~= B.nrow
    error('Wrong dimensions')
end

if strcmpi(A.type,'supermatrix') && strcmpi(B.type,'supermatrix')
    C = hmtx_create('supermatrix',A.irow,B.jcol);
    C.M{1,1} = hmtx_add(hmtx_mult(A.M{1,1},B.M{1,1}),hmtx_mult(A.M{1,2},B.M{2,1}),'+');
    C.M{2,1} = hmtx_add(hmtx_mult(A.M{2,1},B.M{1,1}),hmtx_mult(A.M{2,2},B.M{2,1}),'+');
    C.M{1,2} = hmtx_add(hmtx_mult(A.M{1,1},B.M{1,2}),hmtx_mult(A.M{1,2},B.M{2,2}),'+');
    C.M{2,2} = hmtx_add(hmtx_mult(A.M{2,1},B.M{1,2}),hmtx_mult(A.M{2,2},B.M{2,2}),'+');
    C.eps = max([A.eps B.eps]);
elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'fullmatrix')
    C = hmtx_create('fullmatrix',A.irow,B.jcol);
    C.M = A.M*B.M;
    C.eps = max([A.eps B.eps]);
elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'rkmatrix')
    C = hmtx_create('rkmatrix',A.irow,B.jcol);
    C.U = A.U;
    C.V = B.V*(B.U'*A.V);
    C = rSVD_rkmatrix(C,max([A.eps B.eps]));
    
elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'fullmatrix')
    C = hmtx_create('rkmatrix',A.irow,B.jcol);
    C.U = A.U;
    C.V = B.M'*A.V;
    C = rSVD_rkmatrix(C,max([A.eps B.eps]));
elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'rkmatrix')
    C = hmtx_create('rkmatrix',A.irow,B.jcol);
    C.U = A.M*B.U;
    C.V = B.V;
    C = rSVD_rkmatrix(C,max([A.eps B.eps]));
elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'rkmatrix')
    C = hmtx_create('rkmatrix',A.irow,B.jcol);
    C.U = hmtx_HxM(A,B.U);
    C.V = B.V;
    C = rSVD_rkmatrix(C,max([A.eps B.eps]));
elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'fullmatrix')
    C = hmtx_create('fullmatrix',A.irow,B.jcol);
    C.M = hmtx_HxM(A,B.M);
    C.eps = max([A.eps B.eps]);
elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'supermatrix')
    C = hmtx_create('rkmatrix',A.irow,B.jcol);
    C.V = hmtx_MxH(A.V',B)';
    C.U = A.U;
    C = rSVD_rkmatrix(C,max([A.eps B.eps]));
elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'supermatrix')
    C = hmtx_create('fullmatrix',A.irow,B.jcol);
    C.M = hmtx_MxH(A.M,B);
    C.eps = max([A.eps B.eps]);
end

C = hmtx_compress(C);

end

    
   
    
    