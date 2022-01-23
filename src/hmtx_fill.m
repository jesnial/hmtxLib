function H = hmtx_fill(H,eps,fKern,varargin)

% HMTX_FILL fills a H-matrix created by HMTX_CLUSTER using ACA method for
% far-field blocks
%
% USE:
% H = hmtx_fill(H,eps,fKern,varargin)
%
% INPUTS:
% 'H': H-matrix structure, as created by HMTX_CLUSTER
% 'eps': tolerance for ACA method (use Inf to fill fullmatrix blocks only)
% 'fKern': kernel function passed to ACA
% 'varargin': optional inputs. Available:
%    * 'kMax': maximum rank (default: Inf: unconstrained rank)
%    * 'acaKern': specific kernel for admissible clusters rkmatrix
%      (default: fKern)
%
% OUTPUTS:
% 'H': populated H-matrix
%
% NOTE:
% Recursive function
%
% VERSION:
% Date: 07.01.2013
% Copyright(C) 2013-2020: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 09.01.2013: renamed from FILLHMATRIX
% 14.01.2013: added recompression via SVD
% 15.01.0213: added specific kernel for ACA
% 18.01.2013: added rank to fields
% 22.01.2013: added field eps = 0 to 'fullmatrix' type
% 22.01.2013: added field eps to 'supermatrixmatrix' type
% 16.04.2013: automatic compression of 'rkmatrix' and 'fullmatrix' blocks
% 03.01.2020: added KMAX and ACAKERN as optional inputs

% default values
kMax = Inf;
acaKern = fKern;

% check varargin
for i = 1:2:length(varargin)
    vin = varargin{i};
    if strcmpi(vin,'kMax')
        kMax = varargin{i+1};
    elseif strcmpi(vin,'acaKern')
        acaKern = varargin{i+1};
    else
        warning('MATLAB:hmtx_fill','Wrong VARARGIN parameter, skip ''%s = %s''\n',vin,num2str(varargin{i+1}));
    end
end

% fill matrix
H = fillhmatrix(H,eps,kMax,fKern,acaKern);

% compress 'rkmatrix' and 'fullmatrix' blocks
H = hmtx_compress(H);

end

% recursive function
function H = fillhmatrix(H,eps,kMax,fKern,acakern)

if strcmpi(H.type,'supermatrix')
    H.M{1,1} = fillhmatrix(H.M{1,1},eps,kMax,fKern,acakern);
    H.M{1,2} = fillhmatrix(H.M{1,2},eps,kMax,fKern,acakern);
    H.M{2,1} = fillhmatrix(H.M{2,1},eps,kMax,fKern,acakern);
    H.M{2,2} = fillhmatrix(H.M{2,2},eps,kMax,fKern,acakern);
    H.eps = eps;
elseif strcmpi(H.type,'fullmatrix')
    H.M = fKern(H.irow,H.jcol);
    H.eps = 0;
elseif strcmpi(H.type,'rkmatrix')
    if eps == Inf
        % create empty block
        H.U = sparse([],[],[],H.nrow,0);
        H.V = sparse([],[],[],H.ncol,0);
        H.eps = Inf;
        H.k = 0;
        H.kMax = kMax;
    else
        [H.U,H.V,H.k,H.eps] = aca(acakern,H.irow,H.jcol,eps,kMax);
        H = rSVD_rkmatrix(H,eps,kMax);
    end
end

end
