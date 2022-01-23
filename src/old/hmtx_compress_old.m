function A = hmtx_compress(A,varargin)

% HMTX_COMPRESS loops over the H-matrix and compress 2x2 'rkmatrix' blocks
% and unify 2x2 'fullmatrix' blocks
%
% USE:
% A = hmtx_compress(A)
%
% INPUTS:
% 'A': H-matrix, as created by HMTX_CLUSTER and filled by HMTX_FILL
%
% OUTPUTS:
% 'A': compressed H-matrix
%
% NOTE:
% Recursive function
%
% VERSION:
% Date: 18.01.2013
% Copyright(C) 2013-2020: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 21.01.2013: fixed bug due to wrong reordening
% 02.01.2020: eps set to 0 for fullmatrix

% default
force = 'n';

% check varargin
for i = 1:2:length(varargin)
    vin = varargin{i};
    if strcmpi(vin,'force')
        force = varargin{i+1};
    else
        warning('Wrong VARARGIN parameter, skip ''%s = %s''\n',vin,num2str(varargin{i+1}));
    end
end


if strcmpi(A.type,'supermatrix')
    if strcmpi(A.M{1,1}.type,'supermatrix')
        A.M{1,1} = hmtx_compress(A.M{1,1});
    end
    if strcmpi(A.M{2,1}.type,'supermatrix')
        A.M{2,1} = hmtx_compress(A.M{2,1});
    end
    if strcmpi(A.M{1,2}.type,'supermatrix')
        A.M{1,2} = hmtx_compress(A.M{1,2});
    end
    if strcmpi(A.M{2,2}.type,'supermatrix')
        A.M{2,2} = hmtx_compress(A.M{2,2});
    end
    
    if strcmpi(A.M{1,1}.type,'rkmatrix') && ...
            strcmpi(A.M{1,2}.type,'rkmatrix') && ...
            strcmpi(A.M{2,1}.type,'rkmatrix') && ...
            strcmpi(A.M{2,2}.type,'rkmatrix')
        A = rk_agglomerate(A,force);
    elseif strcmpi(A.M{1,1}.type,'fullmatrix') && ...
            strcmpi(A.M{1,2}.type,'fullmatrix') && ...
            strcmpi(A.M{2,1}.type,'fullmatrix') && ...
            strcmpi(A.M{2,2}.type,'fullmatrix')
        A.type = 'fullmatrix';
        
        % sorting is necessary to be compliant to the upper level matrix indexing
        [A.irow,ix] = sort([A.M{1,1}.irow; A.M{2,1}.irow]);
        [A.jcol,jx] = sort([A.M{1,1}.jcol; A.M{1,2}.jcol]);
        M = [
            A.M{1,1}.M A.M{1,2}.M
            A.M{2,1}.M A.M{2,2}.M
            ];
        A.M = M(ix,jx);
        A.eps = 0;
    end
    
end

end

function A = rk_agglomerate(A,force)

% here A is a 'supermatrix' made by four 'rkmatrix'
% See M. Bebendorf, "Hierarchical matrices", Springer, 2008, pp. 18

if A.eps == Inf
    % create empty block
    A.type = 'rkmatrix';
    A.U = sparse([],[],[],A.nrow,0);
    A.V = sparse([],[],[],A.ncol,0);
    A.eps = Inf;
    A.k = 0;
    A.kMax = [];
    A = rmfield(A,'M');
else
    % QR factorizations
    [Q1,R1] = qr([A.M{1,1}.U A.M{1,2}.U],0);
    [Q2,R2] = qr([A.M{2,1}.U A.M{2,2}.U],0);
    [Q3,R3] = qr([A.M{1,1}.V A.M{2,1}.V],0);
    [Q4,R4] = qr([A.M{1,2}.V A.M{2,2}.V],0);

    % matrix R to be supplied to SVD
    R = [
        R1(:,1:A.M{1,1}.k)*R3(:,1:A.M{1,1}.k)' R1(:,1+A.M{1,1}.k:end)*R4(:,1:A.M{1,2}.k)'
        R2(:,1:A.M{2,1}.k)*R3(:,1+A.M{1,1}.k:end)' R2(:,1+A.M{2,1}.k:end)*R4(:,1+A.M{1,2}.k:end)'
        ];

    % SVD
    [U,S,V] = svd(R,'econ');
    sigma = diag(S);

    % set threshold for truncation
    eps = max([A.M{1,1}.eps, A.M{1,2}.eps, A.M{2,1}.eps, A.M{2,2}.eps]);
    kMax = min([A.M{1,1}.kMax, A.M{1,2}.kMax, A.M{2,1}.kMax, A.M{2,2}.kMax]);
    
    % search rank
    sigma2 = sigma.^2;
    normUSV2 = cumsum(sigma2);
    k = min([kMax find((normUSV2(end)-normUSV2)/normUSV2(end) <= eps^2,1,'first')]);

    % compute storage
    snew = (A.nrow+A.ncol)*k;
    sold = (A.M{1,1}.nrow+A.M{1,1}.ncol)*A.M{1,1}.k+...
        (A.M{2,1}.nrow+A.M{2,1}.ncol)*A.M{2,1}.k+...
        (A.M{1,2}.nrow+A.M{1,1}.ncol)*A.M{1,2}.k+...
        (A.M{2,2}.nrow+A.M{1,1}.ncol)*A.M{2,2}.k;

    % overwrite if agglomeration saves storage
    %if snew < sold || strcmpi(force,'y')

        % sorting is necessary to be compliant with upper level matrices
        [A.irow,ix] = sort([A.M{1,1}.irow; A.M{2,1}.irow]);
        [A.jcol,jx] = sort([A.M{1,1}.jcol; A.M{1,2}.jcol]);

        A.type = 'rkmatrix';
        U = blkdiag(Q1,Q2)*U(:,1:k)*S(1:k,1:k);
        V = blkdiag(Q3,Q4)*V(:,1:k);
        A.U = U(ix,:);
        A.V = V(jx,:);

        A.k = k;
        A.kMax = kMax;
        A.eps = sqrt((normUSV2(end)-normUSV2(k))/normUSV2(end));
        A = rmfield(A,'M');
    %end
end

end

