function Ark = full2rkmatrix(Af,kMax)

% FULL2RKMATRIX converts fullmatrix to rkmatrix via truncated SVD
%
% USE:
% Ark = full2rkmatrix(Af,kMax)
%
% INPUTS:
% 'Af': H-matrix in fullmatrix format
% 'kMax': maximum rank
%
% OUTPUTS:
% 'Ark': H-matrix in rkmatrix format
%
% NOTE:
% See Grasedyck, W. Hackbusch, "Construction and Arithmetics of
% H-matrices", Computing, vol. 70, 2003, p. 305
%
% VERSION:
% Date: 04.02.2014
% Copyright(C) 2014-2020: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 03.01.2020: added KMAX as input

% create structure
Ark = hmtx_create('rkmatrix',Af.irow,Af.jcol);

% SVD
[U,S,V] = svd(Af.M,'econ');
sigma = diag(S);

% search rank
sigma2 = sigma.^2;
normUSV2 = cumsum(sigma2);
r = min([kMax find((normUSV2(end)-normUSV2)/normUSV2(end) <= Af.eps^2,1,'first')]);

% if r == Inf
%     Ark.U = zeros(Af.nrow,0);
%     Ark.V = zeros(Af.ncol,0);
%     Ark.k = 0;
%     Ark.kMax = kMax;
%     Ark.eps = Af.eps;
% else
    
    % load
    Ark.U = U(:,1:r)*S(1:r,1:r);
    Ark.V = V(:,1:r);
    Ark.k = r;
    Ark.kMax = kMax;
    Ark.eps = max([Af.eps,sqrt((normUSV2(end)-normUSV2(r))/normUSV2(end))]);
%end
% % change type
% A.type = 'rkmatrix';
% 
% % SVD
% [U,S,V] = svd(A.M,'econ');
% sigma = diag(S);
% 
% % search rank
% sigma2 = sigma.^2;
% normUSV2 = cumsum(sigma2);
% r = find((normUSV2(end)-normUSV2)/normUSV2(end) <= A.eps^2,1,'first');
% 
% % load
% A.U = U(:,1:r)*S(1:r,1:r);
% A.V = V(:,1:r);
% A.k = r;
% A.eps = max([A.eps,sqrt((normUSV2(end)-normUSV2(r))/normUSV2(end))]);
% 
% % erease the full matrix
% A.M = [];

end
