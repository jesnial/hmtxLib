function As = rkmatrix2super(Ark,As)

% RKMATRIX2SUPER converts rkmatrix to supermatrix according to row/col
% indices
%
% USE:
% As = rkmatrix2super(Ark,irow1,jcol1)
%
% INPUTS:
% 'Ark': H-matrix in rkmatrix format
% 'As': empty H-matrix in supermatrix format
%
% OUTPUTS:
% 'As': H-matrix in supermatrix format
%
% NOTE:
% This function uses a local faster setdiff for sorted vectors
%
% VERSION:
% Date: 24.01.2013
% Copyright(C) 2013-2019: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 05.02.2014: input matrix overwritten
% 07.12.2019: added help
% 07.12.2019: separated input/output
% 10.12.2019: new refurbished routine
% 13.12.2019: added eps field

% identify row index split
[~,irowloc{2}] = setdiff(As.irow,As.M{1,1}.irow);
irowloc{1} = setdiff(1:As.nrow,irowloc{2});

% identify col index split
[~,jcolloc{2}] = setdiff(As.jcol,As.M{1,1}.jcol);
jcolloc{1} = setdiff(1:As.ncol,jcolloc{2});

% load original row/col indices
irow{1} = Ark.irow(irowloc{1});
irow{2} = Ark.irow(irowloc{2});
jcol{1} = Ark.jcol(jcolloc{1});
jcol{2} = Ark.jcol(jcolloc{2});

% populate the supermatrix blocks
for i = 1:2
    for j = 1:2
        % block i,j
        if strcmpi(As.M{i,j}.type,'rkmatrix')
            As.M{i,j} = hmtx_create('rkmatrix',irow{i},jcol{j});
            As.M{i,j}.eps = Ark.eps;
            As.M{i,j}.k = Ark.k;
            As.M{i,j}.U = Ark.U(irowloc{i},:);
            As.M{i,j}.V = Ark.V(jcolloc{j},:);
        elseif strcmpi(As.M{i,j}.type,'supermatrix')
            Atmp = hmtx_create('rkmatrix',irow{i},jcol{j});
            Atmp.eps = Ark.eps;
            Atmp.k = Ark.k;
            Atmp.U = Ark.U(irowloc{i},:);
            Atmp.V = Ark.V(jcolloc{j},:);
            As.M{i,j} = rkmatrix2super(Atmp,As.M{i,j});
        elseif strcmpi(As.M{i,j}.type,'fullmatrix')
            Atmp = hmtx_create('rkmatrix',irow{i},jcol{j});
            Atmp.eps = Ark.eps;
            Atmp.k = Ark.k;
            Atmp.U = Ark.U(irowloc{i},:);
            Atmp.V = Ark.V(jcolloc{j},:);
            As.M{i,j} = rkmatrix2full(Atmp);
        end
    end
end
As.eps = Ark.eps;

% % identify row index split
% [irow2,irowloc2] = setdiff(Ark.irow,irow1);
% irowloc1 = setdiff(1:Ark.nrow,irowloc2);
% % identify col index split
% [jcol2,jcolloc2] = setdiff(Ark.jcol,jcol1);
% jcolloc1 = setdiff(1:Ark.ncol,jcolloc2);
% 
% % create structure
% As = hmtx_create('supermatrix',Ark.irow,Ark.jcol);
% 
% % block 1,1
% As.M{1,1} = hmtx_create('rkmatrix',irow1,jcol1);
% As.M{1,1}.eps = Ark.eps;
% As.M{1,1}.U = Ark.U(irowloc1,:);
% As.M{1,1}.V = Ark.V(jcolloc1,:);
% As.M{1,1}.k = Ark.k;
% 
% % block 1,2
% As.M{1,2} = hmtx_create('rkmatrix',irow1,jcol2);
% As.M{1,2}.eps = Ark.eps;
% As.M{1,2}.U = Ark.U(irowloc1,:);
% As.M{1,2}.V = Ark.V(jcolloc2,:);
% As.M{1,2}.k = Ark.k;
% 
% % block 2,1
% As.M{2,1} = hmtx_create('rkmatrix',irow2,jcol1);
% As.M{2,1}.eps = Ark.eps;
% As.M{2,1}.U = Ark.U(irowloc2,:);
% As.M{2,1}.V = Ark.V(jcolloc1,:);
% As.M{2,1}.k = Ark.k;
% 
% % block 2,2
% As.M{2,2} = hmtx_create('rkmatrix',irow2,jcol2);
% As.M{2,2}.eps = Ark.eps;
% As.M{2,2}.U = Ark.U(irowloc2,:);
% As.M{2,2}.V = Ark.V(jcolloc2,:);
% As.M{2,2}.k = Ark.k;

end

function [c,ic] = setdiff(a,b)

% local fast setdiff for sorted vectors
% hp: a,b with unique sorted numbers

logUA = ~ismembc(a,b);
c = a(logUA);
ic = find(logUA);

end
