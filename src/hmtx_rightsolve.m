function X = hmtx_rightsolve(U,B)

% HMTX_RIGHTSOLVE solves the upper triangular linear problem XU=B
%
% USE:
% X = hmtx_rightsolve(U,B)
%
% INPUTS:
% 'U': upper triangular H-matrix
% 'B': rhs provided as H-matrix
%
% OUTPUTS:
% 'X': solution of XU=B
%
% NOTE:
% X has the same cluster-tree of B
%
% VERSION:
% Date: 18.12.2019
% Copyright(C) 2019-2020: Fabio Freschi (fabio.freschi@polito.it)
%                         Jessie Levillain (jessielevillain@gmail.com)
%
% HISTORY:
%

% preallocate
X = hmtx_copystruct(B);

% call recursive routine

X = rightsolve(U,B,X);
% UT = hmtx_transpose(U);
% BT = hmtx_transpose(B);
% X = hmtx_transpose(hmtx_leftsolve(UT,BT));
end

function X = rightsolve(U,B,X) 

if B.nrow == 0
    B
    B.irow
    B.jcol
elseif strcmpi(U.type,'fullmatrix') && strcmpi(B.type,'fullmatrix')
    % result as fullmatrix
    X = hmtx_create('fullmatrix',B.irow,B.jcol);
    % use slash between full matrices
    X.M = B.M/U.M;
    X.eps = max([U.eps,B.eps]);
elseif strcmpi(U.type,'fullmatrix') && strcmpi(B.type,'rkmatrix')
    % result as rkmatrix
    X = hmtx_create('rkmatrix',B.irow,B.jcol);
    % the slash can be applied to the V matrix only
    X.U = B.U;
    X.V = U.M'\B.V;
    X.k = B.k;
    X.kMax = B.kMax;
    X.eps = max([U.eps,B.eps]);
    % recompress
    %X = rSVD_rkmatrix(X,X.eps,X.kMax);
elseif strcmpi(U.type,'fullmatrix') && strcmpi(B.type,'supermatrix')
    % temporarily convert B to fullmatrix (here it should not be big)
    Bf = hmtx_create('fullmatrix',B.irow,B.jcol);
    Bf = hmtx_changeformat(B,Bf);
    Bf.eps = B.eps;
    % solve as fullmatrix
    Xtmp = hmtx_copystruct(Bf);
    Xtmp = rightsolve(U,Bf,Xtmp);
    % convert back to supermatrix
    X = hmtx_changeformat(Xtmp,X);
elseif strcmpi(U.type,'supermatrix') && ...
        (strcmpi(B.type,'rkmatrix') || strcmpi(B.type,'fullmatrix'))
    % check if B has enough columns to be converted to supermatrix
    if B.nrow > 1
        % convert B to supermatrix
        Bs = hmtx_create('supermatrix',B.irow,B.jcol);
        % load sub-blocks dimentions (any row division is good)
        idx = round(B.nrow/2);
        irow{1} = B.irow(1:idx);
        irow{2} = B.irow(idx+1:end);
        for i = 1:2
            for j = 1:2
                Bs.M{i,j}.type = B.type;
                Bs.M{i,j}.irow = irow{i};
                Bs.M{i,j}.jcol = U.M{i,j}.jcol;
            end
        end
        Bs = hmtx_changeformat(B,Bs);
        % solve
        Xtmp = hmtx_copystruct(Bs);
        Xtmp = rightsolve(U,Bs,Xtmp);
        % convert
        X = hmtx_changeformat(Xtmp,X);
    else
        % B is actually a vector, convert L to full
        Uf = hmtx_create('fullmatrix',U.irow,U.jcol);
        Uf = hmtx_changeformat(U,Uf);
        X = rightsolve(Uf,B,X);
    end
elseif strcmpi(U.type,'supermatrix') && strcmpi(B.type,'supermatrix')
    % solve upper right block to get X11 and X12
    X.M{1,1} = rightsolve(U.M{1,1},B.M{1,1},X.M{1,1});
    X.M{2,1} = rightsolve(U.M{1,1},B.M{2,1},X.M{2,1});
    
    % solve X12*U22 = B12 - X11*U12
    XU = hmtx_copystruct(B.M{1,2});
    % here XU = X11*U12
    XU = hmtx_mult(X.M{1,1},U.M{1,2},XU);
    % here RHS = B12 - X11*U12
    RHS = hmtx_add(B.M{1,2},XU,'-');
    % solve block 2,1 using previous results
    X.M{1,2} = rightsolve(U.M{2,2},RHS,X.M{1,2});
    
    % solve X22*U22 = B22-X21*U12
    XU = hmtx_copystruct(B.M{2,2});
    % here XU = X21*U12
    XU = hmtx_mult(X.M{2,1},U.M{1,2},XU);
    % here RHS = B22-X21*U12
    RHS = hmtx_add(B.M{2,2},XU,'-');
    % solve block 2,2 using previous results
    X.M{2,2} = rightsolve(U.M{2,2},RHS,X.M{2,2});
else
    error('Matlab:hmtx_rightsolve is called with a rkmatrix diagonal block\n');
end

end



% elseif strcmpi(L.type,'supermatrix') && strcmpi(B.type,'fullmatrix')
%     % check on col dimension
%     if B.ncol > 1
%         % convert B to supermatrix
%         Bs = hmtx_create('supermatrix',B.irow,B.jcol);
%         % load sub-blocks dimentions (any col division is good)
%         idx = round(B.ncol/2);
%         jcol{1} = B.jcol(1:idx);
%         jcol{2} = B.jcol(idx+1:end);
%         for i = 1:2
%             for j = 1:2
%                 Bs.M{i,j}.type = B.type;
%                 Bs.M{i,j}.irow = L.M{i,j}.irow;
%                 Bs.M{i,j}.jcol = jcol{j};
%             end
%         end
%         Bs = hmtx_changeformat(B,Bs);
%         % solve
%         Xtmp = hmtx_copystruct(Bs);
%         Xtmp = leftsolve(L,Bs,Xtmp);
%         % convert
%         X = hmtx_changeformat(Xtmp,X);
%     else
%         % convert L to full
%         Lf = hmtx_create('fullmatrix',L.irow,L.jcol);
%         Lf = hmtx_changeformat(L,Lf);
%         X = leftsolve(Lf,B,X);
%     end
% elseif strcmpi(L.type,'rkmatrix')
%     % this should never happen (to be removed)
%     warning('L is rkmatrix\n')
%     Bk = hmtx_create('rkmatrix',B.irow,B.jcol);
%     Bk = hmtx_changeformat(B,Bk);
%     X.U = mldivide(A.U*A.V',Bk.U);
%     X.V =  Bk.V;
% end

% function X = Lsolve(L,B,X)
% 
% if strcmpi(L.type,'fullmatrix')
%     % prepare full rhs
%     Bf = hmtx_create('fullmatrix',B.irow,B.jcol);
%     Bf = hmtx_changeformat(B,Bf);
%     % solve
%     Xtmp = hmtx_create('fullmatrix',L.jcol,B.jcol);
%     Xtmp.M = L.M\Bf.M;
%     Xtmp.eps = max([L.eps,B.eps]);
%     % load on X
%     X = hmtx_changeformat(Xtmp,X);
% elseif strcmpi(L.type,'supermatrix')
%     % convert B to supermatrix
%     Bs = hmtx_copystruct(X);
%     Bs = hmtx_changeformat(B,Bs);
%     
%     % solve upper right block to get X11 and X12
%     X.M{1,1} = Lsolve(L.M{1,1},Bs.M{1,1},X.M{1,1});
%     X.M{1,2} = Lsolve(L.M{1,1},Bs.M{1,2},X.M{1,2});
%     
%     % create matrix structures for multiplications
%     L21X11 = hmtx_copystruct(B.M{2,1});
%     L21X11 = hmtx_mult(L.M{2,1},X.M{1,1},L21X11);
%     RHS1 = hmtx_add(B.M{2,1},L21X11,'-');
%     L21X12 = hmtx_copystruct(B.M{2,2});
%     L21X12 = hmtx_mult(L.M{2,1},X.M{1,2}, L21X12);
%     RHS2 = hmtx_add(B.M{2,2},L21X12,'-');
%     % solve other two blocks using previous results
%     X.M{2,1} = Lsolve(L.M{2,2},RHS1,X.M{2,1});
%     X.M{2,2} = Lsolve(L.M{2,2},RHS2,X.M{2,2});
% elseif strcmpi(L.type,'rkmatrix')  
%     disp('warning\n')
%     Bk = hmtx_create('rkmatrix',B.irow,B.jcol);
%     Bk = hmtx_changeformat(B,Bk);
%     X.U = mldivide(A.U*A.V',Bk.U);
%     X.V =  Bk.V;
% end
  


% if strcmpi(U.type,'fullmatrix')
%     % prepare full rhs
%     Bf = hmtx_create('fullmatrix',B.irow,B.jcol);
%     Bf = formatconversion(B,Bf);
%     % solve
%     Xtmp = hmtx_create('fullmatrix',U.jcol,B.jcol);
%     Xtmp.M = U.M/Bf.M;
%     Xtmp.eps = max([U.eps,B.eps]);
%     % load on X
%     X = formatconversion(Xtmp,X);
% elseif strcmpi(U.type,'supermatrix')
%     % convert B to supermatrix
%     Bs = hmtx_copystruct(U);
%     Bs = formatconversion(B,Bs);
%     % solve upper right block to get X11 and X12
%     X.M{1,1} = Usolve(U.M{1,1},Bs.M{1,1}, X.M{1,1});
%     X.M{2,1} = Usolve(U.M{1,1},Bs.M{2,1}, X.M{2,1});
%     %create matrix structures for multiplications
%     X11U12 = hmtx_copystruct(B.M{1,2});
%     X11U12 = hmtx_mult(X.M{1,1},U.M{1,2}, X11U12);
%     RHS1 = hmtx_add(B.M{1,2}, X11U12, '-');
%     X21U12 = hmtx_copystruct(B.M{2,2});
%     X21U12 = hmtx_mult(X.M{2,1},U.M{1,2}, X21U12);
%     RHS2 = hmtx_add(B.M{2,2}, X21U12, '-');
%     % solve other two blocks using previous results
%     X.M{1,2} = Usolve(U.M{2,2}, RHS1, X.M{1,2});
%     X.M{2,2} = Usolve(U.M{2,2}, RHS2, X.M{2,2});
% elseif strcmpi(U.type,'rkmatrix')  
%     disp('warning')
%     Bk = hmtx_create('rkmatrix',B.irow,B.jcol);
%     Bk = formatconversion(B,Bk);
%     X.U = mldivide(A.U*A.V',Bk.U);
%     X.V =  Bk.V;
% end
% 
% end
% 
% function Hout = formatconversion(Hin,Hout)
% 
% if strcmpi(Hin.type,Hout.type)
%     Hout = Hin;
% else
%     if strcmpi(Hin.type,'supermatrix') && strcmpi(Hout.type,'fullmatrix')
%         Hout = super2full(Hin);
%     elseif strcmpi(Hin.type,'supermatrix') && strcmpi(Hout.type,'rkmatrix')
%         Hout = super2rkmatrix(Hin);
%     elseif strcmpi(Hin.type,'rkmatrix') && strcmpi(Hout.type,'fullmatrix')
%         Hout = rkmatrix2full(Hin);
%     elseif strcmpi(Hin.type,'fullmatrix') && strcmpi(Hout.type,'rkmatrix')
%         Hout = full2rkmatrix(Hin);
%     elseif strcmpi(Hin.type,'rkmatrix') && strcmpi(Hout.type,'supermatrix')
%         Hout = rkmatrix2super(Hin,Hout);
%     elseif strcmpi(Hin.type,'fullmatrix') && strcmpi(Hout.type,'supermatrix')
%         Hout = full2super(Hin,Hout);
%     end
% end
% 
% end
