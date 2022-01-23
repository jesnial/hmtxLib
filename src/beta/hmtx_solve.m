 
function X = hmtx_solve(A,B,side)
% solves matrix system AX=B for X if left
% solves XA = B if right with A triangular
%
% USE:
% X = solve_syst(A,B,side)
%
% INPUTS:
% 'A' :H-matrix upper or lower triangular 
% 'B' : H-matrix
% 'side' : string indicating where 'X' is
% 
%
% OUTPUTS:
% 'X' : H-matrix
%
% NOTE:
% See S. Borm, L. Grasedyck, W. Hackbusch, "Hierarchical Matrices", 2003,
% (revised version: June 2006), pp. 119

%
% VERSION:
% Date: 19.11.2019
% Copyright(C) 2019: Jessie Levillain (jessielevillain@gmail.com)
%
% HISTORY:
%
if strcmpi(side,'right')
    X = hmtx_rightsolve(A,B);
    
elseif strcmpi(side, 'left')
    X = hmtx_leftsolve(A,B);
%     AT =hmtx_transpose(A);
%     BT = hmtx_transpose(B);
%     X = hmtx_transpose(solve_R(AT,BT));
    
end
end 


% function X = solve_L(A,B)
% % solves matrix system LX=B for X 
% %
% % USE:
% % X = solve_L(A,B)
% %
% % INPUTS:
% % 'A' :H-matrix lower triangular
% % 'B' : H-matrix
% % 
% %
% % OUTPUTS:
% % 'X' : H-matrix
% %
% %
% % VERSION:
% % Date: 19.11.2019
% % Copyright(C) 2019: Jessie Levillain (jessielevillain@gmail.com)
% %
% % HISTORY:
% % 26.11.2019 : A and B can be of different types
% %
% if A.nrow ~= B.nrow
%     error('Wrong dimensions')
% end
% 
% if strcmpi(A.type,'supermatrix') &&  strcmpi(B.type,'supermatrix')
%      X = hmtx_create('supermatrix', A.jcol, B.jcol);
%      X.eps = A.eps; %check which precision I should take
%      
%      X.M{1,1} = solve_L(A.M{1,1},B.M{1,1});
%      X.M{1,2} = solve_L(A.M{1,1}, B.M{1,2});
%      L21X11 = hmtx_copystruct(B.M{2,1});
%      L21X11 = hmtx_mult(A.M{2,1}, X.M{1,1}, L21X11);
%      RHS1 = hmtx_add(B.M{2,1}, L21X11, '-');
% %      A.M{2,1}
% %      X.M{1,2}
%     L21X12 = hmtx_copystruct(B.M{2,2});
%     L21X12 = hmtx_mult(A.M{2,1}, X.M{1,2}, L21X12);
%      RHS2 = hmtx_add(B.M{2,2}, L21X12, '-');
%      X.M{2,2} = solve_L(A.M{2,2}, RHS1);
%      X.M{2,1} = solve_L(A.M{2,2},RHS2);
%      
% elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'fullmatrix')
%     X = hmtx_create('fullmatrix',A.jcol, B.jcol);
% %     X.M = A.M\B.M;
% %try using more specific solver
%     X.M = mldivide(A.M, B.M);
%     
% elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'rkmatrix')
%    X = hmtx_create('rkmatrix',A.jcol, B.jcol);
% %   disp('rkrk')
%    % according to hmtx_mult 
%    % B.U = A.U *A.V' * X.U, B.V = X.V
%    %X.U = (A.U*A.V')\B.U;
%    X.U = mldivide(A.U*A.V',B.U);
%    X.V =  B.V;
%    %X
%    %X = rSVD_rkmatrix(X, max(A.eps, B.eps));
%    %X = hmtx_compress(X);
%    
%    
% elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'fullmatrix')
%  %   disp('rkfull')
%     X = hmtx_create('fullmatrix',A.jcol,B.jcol);
%     Af = hmtx_full(A);
%     %X.M = Af.M\B.M;
%     X.M = mldivide(Af,B.M);
%     
% elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'rkmatrix')
%   %  disp('fullrk')
%     X = hmtx_create('rkmatrix',A.jcol,B.jcol);
%     X.U = mldivide(A.M,B.U);
%     %X.U = A.M\B.U;
%     X.V =  B.V;
%     %X = rSVD_rkmatrix(X,max([A.eps B.eps]));
%     %X
%     
% elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'rkmatrix')
% %     B2 = rk2super(B,A.M{1,1}.jcol,[]);
% %     X = solve_L(A,B2);
%     %X = hmtx_compress(X);
%    % disp('superrk')
%     X = hmtx_create('rkmatrix', A.jcol, B.jcol);
%     Af = hmtx_full(A);
%     X.U = mldivide(Af.M, B.U);
%     X.V = B.V;
%     %X
%     
% elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'supermatrix')
% %     A2 = rk2super(A,[],B.M{1,1}.irow);
% %     X = solve_L(A2,B);
%     %X = hmtx_compress(X); 
%   %  disp('rksuper')
%     UA = hmtx_create('fullmatrix', A.irow, A.k);
%     VAT = hmtx_create('fullmatrix', A.k, A.jcol);
%     UA.M = A.U;
%     VAT.M = transpose(A.V);
%     X2 = solve_L(UA, B);
%     X = solve_L(VAT, X2);
%     %X
%     
% elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'fullmatrix')
%     A2 = hmtx_create('fullmatrix', A.irow, A.jcol);
%     A2.M = hmtx_full(A);
%     X = solve_L(A2,B);
%     
% elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'supermatrix')
%     B2 = hmtx_create('fullmatrix', B.irow, B.jcol);
%     B2.M = hmtx_full(B);
%     X = solve_L(A,B2);
% 
% 
% 
% end
% %X = hmtx_compress(X);
% end
% 
% function X = solve_R(A,B)
% % solves matrix system XA=B for X 
% %
% % USE:
% % X = solve_R(A,B)
% %
% % INPUTS:
% % 'A' :H-matrix upper triangular
% % 'B' : H-matrix
% % 
% %
% % OUTPUTS:
% % 'X' : H-matrix
% %
% %
% % VERSION:
% % Date: 19.11.2019
% % Copyright(C) 2019: Jessie Levillain (jessielevillain@gmail.com)
% %
% % HISTORY:
% % 26.11.2019 : A and B can be of different types
% %
% 
% %disp(A.ncol);
% %disp(B.ncol);
% 
% if A.ncol ~= B.ncol
%     error('Wrong dimensions')
%     %disp('warning');
% end
% if strcmpi(A.type,'supermatrix') &&  strcmpi(B.type,'supermatrix')
%      X = hmtx_create('supermatrix', B.irow, A.irow);
%      X.eps = A.eps; %check which precision I should take
%      
%      X.M{1,1} = solve_R(A.M{1,1},B.M{1,1});
%      X.M{2,1} = solve_R(A.M{1,1}, B.M{2,1});
% %      disp('X');
% %      X.M{1,1}
% %      disp('A');
% %      A
% %      A.M{1,2}
% %      A.M{1,1}
% %      A.M{2,2}
% %      A.M{2,1}
%      aaa = hmtx_copystruct(X.M{1,1});
%      aaa = hmtx_mult(X.M{1,1}, A.M{1,2}, aaa);
%      Y = hmtx_add(B.M{1,2},aaa , '-');
%      bbb = hmtx_copystruct(X.M{2,1});
%      bbb = hmtx_mult(X.M{2,1}, A.M{1,2}, bbb)
%      Z = hmtx_add(B.M{2,2},bbb , '-');
%      X.M{1,2} = solve_R(A.M{2,2}, Y);
%      X.M{2,2} = solve_R(A.M{2,2},Z);
%      
% elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'fullmatrix')
%     X = hmtx_create('fullmatrix',B.irow, A.irow);
%     %X.M = A.M/B.M;
%     X.M = mrdivide(A.M,B.M);
%     X.eps = 0;
%     
% elseif strcmpi(A.type,'rkmatrix') && strcmpi(A.type,'rkmatrix')
%    X = hmtx_create('rkmatrix',B.irow, A.irow);
%    % according to hmtx_mult 
%    % B.U = X.U, B.V = A.V*(A.U'*X.V)
%    X.U = B.U;
%    %X.V = (A.V*A.U')/B.V;
%    X.V = mrdivide((A.V*A.U'),B.V);
%    X.eps = max(A.eps, B.eps);
%    %X =  rSVD_rkmatrix(X, max(A.eps, B.eps));
%    %X = hmtx_compress(X)
%   
%    % matrices A and B should have the same structure, check if more elseif
%    % should be added
%   
% elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'fullmatrix')
%     X = hmtx_create('fullmatrix',B.irow,A.irow);
%     Af = hmtx_full(A);
%     X.M = mrdivide(Af,B.M);
%     X.eps = A.eps;
% %     X.M = Af/B.M;
%     
% elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'rkmatrix')
%     X = hmtx_create('rkmatrix',B.irow,A.irow);
%     X.U = B.U;
% %     A.irow
% %     A.jcol
% %     size(B.V)
% %     A.M
% %     B.V
% %     B.jcol
% %     B.irow
% %    X.V =  A.M\B.V;
%     X.V = mldivide (A.M, B.V);
%     X.eps = B.eps;
% %     X
% %     X = rSVD_rkmatrix(X,max([A.eps B.eps]));
%     
% elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'rkmatrix')
%     %B2 = rk2super(B,A.M{1,1}.jcol,[]);
% %     B2 = rk2super(B,[],[]);
% %     X = solve_R(A,B2);
%     %X = hmtx_compress(X);
%     
%     X = hmtx_create('rkmatrix', B.irow, A.irow);
%     X.U = B.U;
%     AT = hmtx_full(hmtx_transpose(A));
%     X.V = mldivide(AT.M, B.V);
%     X.eps = B.eps;
% %     X = rSVD_rkmatrix(X,max([A.eps B.eps]));
%     
% elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'supermatrix')
% %     A2 = rk2super(A,[],B.M{1,1}.irow);
% %     X = solve_R(A2,B);
%     %X = hmtx_compress(X);
%     UA = hmtx_create('fullmatrix', A.irow, A.k);
%     VAT = hmtx_create('fullmatrix', A.k, A.jcol);
%     UA.M = A.U;
%     VAT.M = transpose(A.V);
%     X2 = solve_R(VAT, B);
%     X = solve_R(UA, X2);
%     X.eps = A.eps;
%     
% elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'fullmatrix')
%     %B2 = full2super(B,A.M{1,1}.jcol,[]);
%     A2 = hmtx_create('fullmatrix', A.irow, A.jcol)
%     A2.M = hmtx_full(A);
%     X = solve_R(A2,B);
%     X.eps = max(A.eps, B.eps);
%     %X = hmtx_compress(X);
% 
%     
% elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'supermatrix')
%     %A2 = full2super(A,[],B.M{1,1}.irow);
%     B2 = hmtx_create('fullmatrix', B.irow, B.jcol)
%     B2.M = hmtx_full(B);
%     X = solve_R(A,B2);
%     X.eps = 0;
%     %X = hmtx_compress(X);
% 
% 
% end
% %X = hmtx_compress(X);
% end