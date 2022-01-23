  function [Lh,Uh] = hmtx_lu(A)

% Calculates an (approximate) H-LU decomposition of A such that
% A ~ Lh*Uh
%
% USE:
% [Lh,Uh, k, err] = hmtx_lu(A)
%
% INPUTS:
% 'A' :H-matrix
% 
%
% OUTPUTS:
% 'Lh': matrix with rows
% 'Uh': matrix with columns
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

 
    
if strcmpi(A.type,'supermatrix')
     Lh = hmtx_create('supermatrix', A.irow, A.jcol);
     Uh = hmtx_create('supermatrix', A.irow, A.jcol);
     Lh.eps = A.eps;
     Uh.eps = A.eps;
     
     %create blocks full of 0 to avoid problems
     Lh.M{1,2} = hmtx_create('rkmatrix', A.M{1,2}.irow, A.M{1,2}.jcol);
     Lh.M{1,2}.eps = 0;
     Lh.M{1,2}.k = 1;
     Lh.M{1,2}.U = zeros(A.M{1,2}.nrow, Lh.M{1,2}.k);
     Lh.M{1,2}.V = zeros(A.M{1,2}.ncol, Lh.M{1,2}.k);
     
     Uh.M{2,1} = hmtx_create('rkmatrix', A.M{2,1}.irow, A.M{2,1}.jcol);
     Uh.M{2,1}.eps = 0;
     Uh.M{2,1}.k = 1;
     Uh.M{2,1}.U = zeros(A.M{2,1}.nrow, Uh.M{2,1}.k);
     Uh.M{2,1}.V = zeros(A.M{2,1}.ncol, Uh.M{2,1}.k);
     
     %find LU for first block
     [Lh.M{1,1},Uh.M{1,1}] = hmtx_lu(A.M{1,1});
     
     %solve for diagonal blocs
     Uh.M{1,2} = hmtx_solve(Lh.M{1,1}, A.M{1,2} , 'left');
     Lh.M{2,1} = hmtx_solve(Uh.M{1,1}, A.M{2,1}, 'right');
     
     % find LU for last blocks 
     % L22U22 = A22 - L21U12
     CCC = hmtx_copystruct(A.M{2,2});
     CCC = hmtx_mult(Lh.M{2,1},Uh.M{1,2}, CCC);
     [Lh.M{2,2},Uh.M{2,2}] = hmtx_lu(hmtx_add(A.M{2,2},CCC,'-') ); 

     
elseif strcmpi(A.type,'fullmatrix')
    %use standard LU
    Lh = hmtx_create('fullmatrix',A.irow,A.jcol);
    Uh = hmtx_create('fullmatrix',A.irow,A.jcol);
   [Lh.M,Uh.M] = lu(A.M);
 
% elseif strcmpi(A.type,'rkmatrix')
%    Lh = hmtx_create('rkmatrix',A.irow,A.jcol);
%    Uh = hmtx_create('rkmatrix',A.irow,A.jcol);
%    Lh.eps = A.eps;
%    Uh.eps = A.eps;
%    
%    %find LU for ACA of A
%    [Lh,Uh] = lu(A.U*A.V'); % any better options ?
%    
%    %ACA method and recompression for L and U
%    %[A.Lh.U,A.Lh.V, Lh.k, Lh.eps]=aca(Lh.fkern, Lh.irow, Lh.jcol, Lh.eps); 
%    Lh = rSVD_rkmatrix(A.Lh,A.Lh.eps);
%    Lh = hmtx_compress(Lh); %needed?
%    
%    %[A.Uh.U,A.Uh.V, Uh.k, Uh.eps]=aca(Uh.fkern, Uh.irow, Uh.jcol, Uh.eps);
%    A.Uh = rSVD_rkmatrix(A.Uh,A.Uh.eps);
%    Uh = hmtx_compress(Uh);
   
end

% Lh = hmtx_compress(Lh);
% Uh = hmtx_compress(Uh);
end
