function C = hmtx_mult(A,B,C)

% HMTX_MULT multiplies two H-matrices with the same cluster tree
%
% USE:
% C = hmtx_mult(A,B)
%
% INPUTS:
% 'A': H-matrix, as created by HMTX_CLUSTER and filled by HMTX_FILL
% 'B': H-matrix, as created by HMTX_CLUSTER and filled by HMTX_FILL
%
% OUTPUTS:
% 'C': H-matrix such that C = A*B
%
% NOTE:
% When eps == Inf, the resulting matrix suffers the fill-in of rkmatrices
% blocks due to out-of-diagonal fullmatrix blocks.
%
% VERSION:
% Date: 14.01.2013
% Copyright(C) 2013-2020: Fabio Freschi (fabio.freschi@polito.it)
%                         Jessie Levillain (jessielevillain@gmail.com)
%
% HISTORY:
% 23.01.2013: uses HMTX_CREATE to create H-matrix structure
% 29.11.2019: fixed cross-type cases. supermatrix and rk becomes rk
% 03.12.2019: added recompression
% 16.12.2019: new refurbished version
% 02.01.2020: added k and eps for rkmatrices

if strcmpi(C.type,'supermatrix')
    if strcmpi(A.type,'supermatrix') && strcmpi(B.type,'supermatrix')
        for i = 1:2
            for j = 1:2
                Ca = hmtx_copystruct(C.M{i,j});
                Ca = hmtx_mult(A.M{i,1},B.M{1,j},Ca);
                Cb = hmtx_copystruct(C.M{i,j});
                Cb = hmtx_mult(A.M{i,2},B.M{2,j},Cb);
                C.M{i,j} = hmtx_add(Ca,Cb,'+');
            end
        end
        C.eps = max([A.eps B.eps]);
    elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'supermatrix')
        % solution as fullmatrix
        Ctmp = hmtx_create('fullmatrix',A.irow,B.jcol);
        Ctmp.M = hmtx_MxH(A.M,B);
        Ctmp.eps = max([A.eps B.eps]);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'supermatrix')
        % solution as rkmatrix
        Ctmp = hmtx_create('rkmatrix',A.irow,B.jcol);
        Ctmp.V = hmtx_MxH(A.V',B)';
        Ctmp.U = A.U;
        Ctmp.k = A.k;
        Ctmp.kMax = A.kMax;
        Ctmp.eps = max([A.eps B.eps]);
        Ctmp = rSVD_rkmatrix(Ctmp,Ctmp.eps,Ctmp.kMax);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'fullmatrix')
        % solution as fullmatrix
        Ctmp = hmtx_create('fullmatrix',A.irow,B.jcol);
        Ctmp.M = hmtx_HxM(A,B.M);
        Ctmp.eps = max([A.eps B.eps]);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'fullmatrix')
        % solution as fullmatrix
        Ctmp = hmtx_create('fullmatrix',A.irow,B.jcol);
        Ctmp.M = A.M*B.M;
        Ctmp.eps = max([A.eps B.eps]);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'fullmatrix')
        % solution as rkmatrix
        Ctmp = hmtx_create('rkmatrix',A.irow,B.jcol);
        Ctmp.U = A.U;
        Ctmp.V = B.M'*A.V;
        Ctmp.k = A.k;
        Ctmp.kMax = A.kMax;
        Ctmp.eps = max([A.eps B.eps]);
        Ctmp = rSVD_rkmatrix(Ctmp,Ctmp.eps,Ctmp.kMax);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'rkmatrix')
        % solution as rkmatrix
        Ctmp = hmtx_create('rkmatrix',A.irow,B.jcol);
        Ctmp.U = hmtx_HxM(A,B.U);
        Ctmp.V = B.V;
        Ctmp.k = B.k;
        Ctmp.kMax = B.kMax;
        Ctmp.eps = max([A.eps B.eps]);
        Ctmp = rSVD_rkmatrix(Ctmp,Ctmp.eps,Ctmp.kMax);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'rkmatrix')
        % solution as rkmatrix
        Ctmp = hmtx_create('rkmatrix',A.irow,B.jcol);
        Ctmp.U = A.M*B.U;
        Ctmp.V = B.V;
        Ctmp.k = B.k;
        Ctmp.kMax = B.kMax;
        Ctmp.eps = max([A.eps B.eps]);
        Ctmp = rSVD_rkmatrix(Ctmp,Ctmp.eps,Ctmp.kMax);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'rkmatrix')
        % solution as rkmatrix
        Ctmp = hmtx_create('rkmatrix',A.irow,B.jcol);
        Ctmp.U = A.U;
        Ctmp.V = B.V*(B.U'*A.V);
        Ctmp.k = A.k;
        Ctmp.kMax = min([A.kMax B.kMax]);
        Ctmp.eps = max([A.eps B.eps]);
        Ctmp = rSVD_rkmatrix(Ctmp,Ctmp.eps,Ctmp.kMax);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    end
else
    if strcmpi(A.type,'supermatrix') && strcmpi(B.type,'supermatrix')
        Ctmp = hmtx_create('supermatrix',A.irow,B.jcol);
        Ctmp.eps = max([A.eps B.eps]);
        for i = 1:2
            for j = 1:2
                Ca = hmtx_create(C.type,A.irow,B.jcol);
                Ca = hmtx_mult(A.M{i,1},B.M{1,j},Ca);
                
                Cb = hmtx_create(C.type,A.irow,B.jcol);
                Cb = hmtx_mult(A.M{i,2},B.M{2,j},Cb);
                Ctmp.M{i,j} = hmtx_add(Ca,Cb,'+');
            end
        end
        Ctmp.eps = max([Ctmp.M{1,1}.eps,Ctmp.M{1,2}.eps,Ctmp.M{2,1}.eps,Ctmp.M{2,2}.eps]);
        C = hmtx_changeformat(Ctmp,C);
        % C.eps = max([A.eps B.eps]);
    elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'supermatrix')
        % solution as fullmatrix
        Ctmp = hmtx_create('fullmatrix',A.irow,B.jcol);
        Ctmp.M = hmtx_MxH(A.M,B);
        Ctmp.eps = max([A.eps B.eps]);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'supermatrix')
        % solution as rkmatrix
        Ctmp = hmtx_create('rkmatrix',A.irow,B.jcol);
        Ctmp.V = hmtx_MxH(A.V',B)';
        Ctmp.U = A.U;
        Ctmp.k = A.k;
        Ctmp.kMax = A.kMax;
        Ctmp.eps = max([A.eps B.eps]);
        Ctmp = rSVD_rkmatrix(Ctmp,Ctmp.eps,Ctmp.kMax);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'fullmatrix')
        % solution as fullmatrix
        Ctmp = hmtx_create('fullmatrix',A.irow,B.jcol);
        Ctmp.M = hmtx_HxM(A,B.M);
        Ctmp.eps = max([A.eps B.eps]);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'fullmatrix')
        % solution as fullmatrix
        Ctmp = hmtx_create('fullmatrix',A.irow,B.jcol);
        Ctmp.M = A.M*B.M;
        Ctmp.eps = max([A.eps B.eps]);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'fullmatrix')
        % solution as rkmatrix
        Ctmp = hmtx_create('rkmatrix',A.irow,B.jcol);
        Ctmp.U = A.U;
        Ctmp.V = B.M'*A.V;
        Ctmp.k = A.k;
        Ctmp.kMax = A.kMax;
        Ctmp.eps = max([A.eps B.eps]);
        Ctmp = rSVD_rkmatrix(Ctmp,Ctmp.eps,Ctmp.kMax);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'rkmatrix')
        % solution as rkmatrix
        Ctmp = hmtx_create('rkmatrix',A.irow,B.jcol);
        Ctmp.U = hmtx_HxM(A,B.U);
        Ctmp.V = B.V;
        Ctmp.k = B.k;
        Ctmp.kMax = B.kMax;
        Ctmp.eps = max([A.eps B.eps]);
        Ctmp = rSVD_rkmatrix(Ctmp,Ctmp.eps,Ctmp.kMax);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'rkmatrix')
        % solution as rkmatrix
        Ctmp = hmtx_create('rkmatrix',A.irow,B.jcol);
        Ctmp.U = A.M*B.U;
        Ctmp.V = B.V;
        Ctmp.k = B.k;
        Ctmp.kMax = B.kMax;
        Ctmp.eps = max([A.eps B.eps]);
        Ctmp = rSVD_rkmatrix(Ctmp,Ctmp.eps,Ctmp.kMax);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'rkmatrix')
        % solution as rkmatrix
        Ctmp = hmtx_create('rkmatrix',A.irow,B.jcol);
        Ctmp.U = A.U;
        Ctmp.V = B.V*(B.U'*A.V);
        Ctmp.k = A.k;
        Ctmp.kMax = min([A.kMax B.kMax]);
        Ctmp.eps = max([A.eps B.eps]);
        Ctmp = rSVD_rkmatrix(Ctmp,Ctmp.eps,Ctmp.kMax);
        % conversion to supermatrix
        C = hmtx_changeformat(Ctmp,C);
    end
end

% C = hmtx_compress(C);

end