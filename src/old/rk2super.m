%function A2 = rk2super(A,irow1,jcol1)
function A2 = rk2super(A,nrow1,ncol1)

% 24.01.2013

A2 = hmtx_create('supermatrix',A.irow,A.jcol);
A2.eps = A.eps;

if ~isempty(nrow1)
    irow1 = linspace(1, nrow1, nrow1);
%     [irow2,irowloc2] = setdiff(A.irow,irow1);
%     irowloc1 = setdiff(1:A.nrow,irowloc2);
    irow1 = A.irow(1:nrow1);
    irow2 = A.irow(nrow1+1:A.nrow);
    irowloc1 = 1:nrow1;
    irowloc2 = nrow1+1:A.nrow;
else
    ncut = round(A.nrow/2);
    irow1 = A.irow(1:ncut);
    irow2 = A.irow(ncut+1:end);
    irowloc1 = 1:ncut;
    irowloc2 = ncut+1:A.nrow;    
end
if ~isempty(ncol1)
    jcol1 = linspace(1, ncol1, ncol1);
%     [jcol2,jcolloc2] = setdiff(A.jcol,jcol1);
%     jcolloc1 = setdiff(1:A.ncol,jcolloc2);
    jcol1 = A.jcol(1:ncol1);
    jcol2 = A.jcol(ncol1+1:A.ncol);
    jcolloc1 = 1:ncol1;
    jcolloc2 = ncol1+1:A.ncol;
else 
    ncut = round(A.ncol/2);
    jcol1 = A.jcol(1:ncut);
    jcol2 = A.jcol(ncut+1:end);
    jcolloc1 = 1:ncut;
    jcolloc2 = ncut+1:A.ncol;    
end

% jcol1
% jcol2
% irow1
% irow2

% block 1,1
A2.M{1,1} = hmtx_create('rkmatrix',irow1,jcol1);
A2.M{1,1}.eps = A.eps;
A2.M{1,1}.U = A.U(irowloc1,:);
A2.M{1,1}.V = A.V(jcolloc1,:);
A2.M{1,1}.k = A.k;
% A2.M{1,1}
% A.V(jcolloc1,:)
% A.V
% A

% block 1,2
A2.M{1,2} = hmtx_create('rkmatrix',irow1,jcol2);
A2.M{1,2}.eps = A.eps;
A2.M{1,2}.U = A.U(irowloc1,:);
A2.M{1,2}.V = A.V(jcolloc2,:);
A2.M{1,2}.k = A.k;
%A.M{1,2}

% block 2,1
A2.M{2,1} = hmtx_create('rkmatrix',irow2,jcol1);
A2.M{2,1}.eps = A.eps;
A2.M{2,1}.U = A.U(irowloc2,:);
A2.M{2,1}.V = A.V(jcolloc1,:);
A2.M{2,1}.k = A.k;
%A2.M{2,1}

% block 2,2
A2.M{2,2} = hmtx_create('rkmatrix',irow2,jcol2);
A2.M{2,2}.eps = A.eps;
A2.M{2,2}.U = A.U(irowloc2,:);
A2.M{2,2}.V = A.V(jcolloc2,:);
A2.M{2,2}.k = A.k;
%A2.M{2,2}

end
