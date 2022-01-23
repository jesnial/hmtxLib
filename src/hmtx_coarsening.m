function A = hmtx_coarsening(A,eps,varargin)
    
% HMTX_COARSENING coarsens the H-matrix to a new tolerance or new maximum
% rank
%
% USE:
% A = hmtx_coarsening(A,eps,varargin)
%
% INPUTS:
% 'A': H-matrix structure, as created by HMTX_CLUSTER and filled by HMTX_FILL
% 'eps': tolerance for the Frobenius norm
% 'varargin': optional inputs. Available:
%    * 'kMax': maximum rank (default: Inf: unconstrained rank)
%
% OUTPUTS:
% 'A': coarsened H-matrix
%
% NOTE:
% See S. Borm, L. Grasedyck, W. Hackbusch, "Hierarchical Matrices", 2003,
% (revised version: June 2006), pp. 69
%
% VERSION:
% Date: 24.01.2013
% Copyright(C) 2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 08.01.2020: added VARARGIN

% default values
kMax = Inf;

% check varargin
for i = 1:2:length(varargin)
    vin = varargin{i};
    if strcmpi(vin,'kMax')
        kMax = varargin{i+1};
    else
        warning('MATLAB:hmtx_fill','Wrong VARARGIN parameter, skip ''%s = %s''\n',vin,num2str(varargin{i+1}));
    end
end

if eps > A.eps || eps == Inf
    % call effective routine
    A = coarsening(A,eps,kMax);

    % compress resulting matrix
    A = hmtx_compress(A);
    A.eps = eps;
else
    fprintf('eps > A.eps: no action required\n');
end
    
end

function A = coarsening(A,eps,kMax)
    
if strcmpi(A.type,'supermatrix')
    % recursive call
    A.M{1,1} = coarsening(A.M{1,1},eps,kMax);
    A.M{1,2} = coarsening(A.M{1,2},eps,kMax);
    A.M{2,1} = coarsening(A.M{2,1},eps,kMax);
    A.M{2,2} = coarsening(A.M{2,2},eps,kMax);
elseif strcmpi(A.type,'rkmatrix')
    if eps == Inf
        % delete the block
        A.U = sparse([],[],[],A.nrow,1);
        A.V = sparse([],[],[],A.ncol,1);
        A.eps = Inf;
        A.k = 0;
    else
        % norm values from r = 1 to r = k
        normu2v2 = (sum(A.U.^2,1).*sum(A.V.^2,1));

        % error on norm
        err = normu2v2./cumsum(normu2v2);
        k = min([kMax,find(err <= eps^2,1,'first')]);

        % overwrite if rank is reduced
        if k < A.k
            A.k = k;
            A.U = A.U(:,1:k);
            A.V = A.V(:,1:k);
            A.eps = sqrt(err(k));
        end
        % in any case, set kMax
        A.kMax = kMax;
    end
end
    
end
