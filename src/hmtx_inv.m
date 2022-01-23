function Hinv = hmtx_inv(H)

% HMTX_INV calculate the H-matrix inverse with the same cluster tree of the
% original
%
% USE:
% Hinv = hmtx_inv(H)
%
% INPUTS:
% 'H': H-matrix
%
% OUTPUTS:
% 'Hinv': inverse H-matrix
%
% NOTE:
% See S. Borm, L. Grasedyck, W. Hackbusch, "Hierarchical Matrices", 2003,
% (revised version: June 2006), pp. 87
%
% VERSION:
% Date: 22.01.2013
% Copyright(C) 2013-2020: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 23.01.2013: uses HMTX_CREATE to create H-matrix structure
% 03.12.2020: fixed bug on dimensions of target matrices 
% 05.01.2020: removed unused cases

if strcmpi(H.type,'supermatrix')
    
    % structure
    Hinv = hmtx_create('supermatrix',H.jcol,H.irow);
    Hinv.eps = H.eps;

    % upper left diagonal block inversion
    H11inv = hmtx_inv(H.M{1,1});
    
    % A = inv(H11)*H12
    A = hmtx_copystruct(H.M{1,2});
    A = hmtx_mult(H11inv,H.M{1,2},A);
    
    % S = H22-H21*A
    H21A = hmtx_copystruct(H.M{2,2});
    H21A = hmtx_mult(H.M{2,1},A,H21A);
    S = hmtx_add(H.M{2,2},H21A,'-');
    % inverse of S
    Sinv = hmtx_inv(S);
       
    % B = H21*inv(H11)
    B = hmtx_copystruct(H.M{2,1});
    B = hmtx_mult(H.M{2,1},H11inv,B);
    
    % C = inv(S)*B
    C = hmtx_copystruct(B);
    C = hmtx_mult(Sinv,B,C);
    
    % fill the blocks
    AC = hmtx_copystruct(H11inv);
    AC = hmtx_mult(A,C,AC);
    Hinv.M{1,1} = hmtx_add(H11inv,AC,'+');
    ASinv = hmtx_copystruct(A);
    ASinv = hmtx_mult(A,Sinv,ASinv);
    Hinv.M{1,2} = hmtx_tH(-1,ASinv);
    Hinv.M{2,1} = hmtx_tH(-1,C);
    Hinv.M{2,2} = Sinv;
        
elseif strcmpi(H.type,'fullmatrix')
    Hinv = hmtx_create('fullmatrix',H.jcol,H.irow);
    Hinv.M = inv(H.M);
    Hinv.eps = H.eps;
else
    error('Matlab:hmtx_inv is called with a rkmatrix diagonal block\n');
end

end

