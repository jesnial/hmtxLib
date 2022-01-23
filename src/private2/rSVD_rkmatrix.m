function A = rSVD_rkmatrix(A,eps,kMax)

% RSVD_RKMATRIX perform recompression of a RKMATRIX via truncated SVD
%
% USE:
% A = rsvd_rkmatrix(A,eps)
%
% INPUTS:
% 'A': rkmatrix
% 'eps': tolerance used for ACA method
%
% OUTPUTS:
% 'A': recompressed rkmatrix
%
% NOTE:
% See S. Borm, L. Grasedyck, W. Hackbusch, "Hierarchical Matrices", 2003,
% (revised version: June 2006), pp. 108
%
% VERSION:
% Date: 13.01.2013
% Copyright(C) 2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 14.01.2013: economy size QR and SVD decomposition
% 22.01.2013: correction of EPS calculation
% 30.12.2019: added possibility of A.eps = Inf

if eps ~= Inf
    % recompression with SVD
    [QU,RU] = qr(A.U,0);
    [QV,RV] = qr(A.V,0);
    [U,S,V] = svd(RU*RV','econ');
    sigma = diag(S);
    
    % search rank
    sigma2 = sigma.^2;
    normUSV2 = cumsum(sigma2);
    
    if all(normUSV2 == 0) || isempty(normUSV2)
        A.U = zeros(A.nrow,0);
        A.V = zeros(A.ncol,0);
        A.k = 0;
        A.eps = 0;
        A.kMax = [];
        warning('rSVD_rkmatrix: normUSV2 == 0\n');
    else
        % rank
        k = min([kMax,find((normUSV2(end)-normUSV2)/normUSV2(end) <= eps^2,1,'first')]);
        % overwrite
        A.U = QU*U(:,1:k)*S(1:k,1:k);
        A.V = QV*V(:,1:k);
        A.k = k;
        A.kMax = kMax;
        A.eps = max([eps,sqrt((normUSV2(end)-normUSV2(k))/normUSV2(end))]);
    end
end

end