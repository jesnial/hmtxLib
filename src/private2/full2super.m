function As = full2super(Af,As)

% FULL2SUPER converts fullmatrix to supermatrix according to row/col
% indices
%
% USE:
% As = full2super(Af,As)
%
% INPUTS:
% 'Af': H-matrix in fullmatrix format
% 'As': empty H-matrix in supermatrix format
%
% OUTPUTS:
% 'As': H-matrix in supermatrix format
%
% NOTE:
% This function uses a local faster setdiff for sorted vectors
%
% VERSION:
% Date: 23.01.2013
% Copyright(C) 2013-2019: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 05.02.2014: input matrix overwritten
% 07.12.2019: removed the possibility of empty inputs
% 10.12.2019: new refurbished routine

% identify row index split
[~,irowloc{2}] = setdiff(As.irow,As.M{1,1}.irow);
irowloc{1} = setdiff(1:As.nrow,irowloc{2});

% identify col index split
[~,jcolloc{2}] = setdiff(As.jcol,As.M{1,1}.jcol);
jcolloc{1} = setdiff(1:As.ncol,jcolloc{2});

% load original row/col indices
irow{1} = Af.irow(irowloc{1});
irow{2} = Af.irow(irowloc{2});
jcol{1} = Af.jcol(jcolloc{1});
jcol{2} = Af.jcol(jcolloc{2});

% populate the supermatrix blocks
for i = 1:2
    for j = 1:2
        % block i,j
        if strcmpi(As.M{i,j}.type,'fullmatrix')
            As.M{i,j} = hmtx_create('fullmatrix',irow{i},jcol{j});
            As.M{i,j}.eps = Af.eps;
            As.M{i,j}.M = Af.M(irowloc{i},jcolloc{j});
        elseif strcmpi(As.M{i,j}.type,'supermatrix')
            Atmp = hmtx_create('fullmatrix',irow{i},jcol{j});
            Atmp.M = Af.M(irowloc{i},jcolloc{j});
            Atmp.eps = 0;
            As.M{i,j} = full2super(Atmp,As.M{i,j});
        elseif strcmpi(As.M{i,j}.type,'rkmatrix')
            Atmp = hmtx_create('fullmatrix',irow{i},jcol{j});
            Atmp.M = Af.M(irowloc{i},jcolloc{j});
            Atmp.eps = 0;
            As.M{i,j} = full2rkmatrix(Atmp,As.M{i,j}.kMax);
        end
    end
end

end

function [c,ic] = setdiff(a,b)

% local fast setdiff for sorted vectors
% hp: a,b with unique sorted numbers

logUA = ~ismembc(a,b);
c = a(logUA);
ic = find(logUA);

end

