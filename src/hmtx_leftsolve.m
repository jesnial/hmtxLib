function X = hmtx_leftsolve(L,B)

% HMTX_LEFTSOLVE solves the lower triangular linear problem LX=B
%
% USE:
% X = hmtx_leftsolve(L,B)
%
% INPUTS:
% 'L': lower triangular H-matrix
% 'B': rhs provided as H-matrix
%
% OUTPUTS:
% 'X': solution of LX=B
%
% NOTE:
% X has the same cluster-tree of B
%
% VERSION:
% Date: 17.12.2019
% Copyright(C) 2019-2020: Fabio Freschi (fabio.freschi@polito.it)
%                         Jessie Levillain (jessielevillain@gmail.com)
%
% HISTORY:
% 04.01.2020: refurbished routine
% 05.01.2020: removed unused cases

% solution has the same structure of the RHS
X = hmtx_copystruct(B);

% call recursive routine
X = leftsolve(L,B,X);

end

function X = leftsolve(L,B,X)

if strcmpi(L.type,'fullmatrix') && strcmpi(B.type,'fullmatrix')
    % result as fullmatrix
    X = hmtx_create('fullmatrix',B.irow,B.jcol);
    % use backslash between full matrices
    X.M = L.M\B.M;
    X.eps = max([L.eps,B.eps]);
elseif strcmpi(L.type,'fullmatrix') && strcmpi(B.type,'rkmatrix')
    % result as rkmatrix
    X = hmtx_create('rkmatrix',B.irow,B.jcol);
    % the backslash can be applied to the U matrix only
    X.U = L.M\B.U;
    X.V = B.V;
    X.k = B.k;
    X.kMax = B.kMax;
    X.eps = max([L.eps,B.eps]);
    % recompress
    %X = rSVD_rkmatrix(X,X.eps,X.kMax);
elseif strcmpi(L.type,'fullmatrix') && strcmpi(B.type,'supermatrix')
    % temporarily convert B to fullmatrix (here it should not be big)
    Bf = hmtx_create('fullmatrix',B.irow,B.jcol);
    Bf = hmtx_changeformat(B,Bf);
    Bf.eps = B.eps;
    % solve as fullmatrix
    Xtmp = hmtx_copystruct(Bf);
    Xtmp = leftsolve(L,Bf,Xtmp);
    % convert back to supermatrix
    X = hmtx_changeformat(Xtmp,X);
elseif strcmpi(L.type,'supermatrix') && ...
        (strcmpi(B.type,'rkmatrix') || strcmpi(B.type,'fullmatrix'))
    % check if B has enough columns to be converted to supermatrix
    if B.ncol > 1
        % convert B to supermatrix
        Bs = hmtx_create('supermatrix',B.irow,B.jcol);
        % load sub-blocks dimentions (any col division is good)
        idx = round(B.ncol/2);
        jcol{1} = B.jcol(1:idx);
        jcol{2} = B.jcol(idx+1:end);
        for i = 1:2
            for j = 1:2
                Bs.M{i,j}.type = B.type;
                Bs.M{i,j}.irow = L.M{i,j}.irow;
                Bs.M{i,j}.jcol = jcol{j};
            end
        end
        Bs = hmtx_changeformat(B,Bs);
        % solve
        Xtmp = hmtx_copystruct(Bs);
        Xtmp = leftsolve(L,Bs,Xtmp);
        % convert
        X = hmtx_changeformat(Xtmp,X);
    else
        % B is actually a vector, convert L to full
        Lf = hmtx_create('fullmatrix',L.irow,L.jcol);
        Lf = hmtx_changeformat(L,Lf);
        X = leftsolve(Lf,B,X);
    end
elseif strcmpi(L.type,'supermatrix') && strcmpi(B.type,'supermatrix')
    % solve upper right block to get X11 and X12
    X.M{1,1} = leftsolve(L.M{1,1},B.M{1,1},X.M{1,1});
    X.M{1,2} = leftsolve(L.M{1,1},B.M{1,2},X.M{1,2});
    
    % solve L22*X21 = B21-L21*X11
    LX = hmtx_copystruct(B.M{2,1});
    % here LX = L21*X11
    LX = hmtx_mult(L.M{2,1},X.M{1,1},LX);
    % here RHS = B21-L21*X11
    RHS = hmtx_add(B.M{2,1},LX,'-');
    % solve block 2,1 using previous results
    X.M{2,1} = leftsolve(L.M{2,2},RHS,X.M{2,1});
    
    % solve L22*X22 = B22-L21*X12
    LX = hmtx_copystruct(B.M{2,2});
    % here LX = L21*X12
    LX = hmtx_mult(L.M{2,1},X.M{1,2},LX);
    % here RHS = B22-L21*X12
    RHS = hmtx_add(B.M{2,2},LX,'-');
    % solve block 2,2 using previous results
    X.M{2,2} = leftsolve(L.M{2,2},RHS,X.M{2,2});
else
    error('Matlab:hmtx_leftsolve is called with a rkmatrix diagonal block\n');
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
