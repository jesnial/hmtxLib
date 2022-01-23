function y = hmtx_mvm(H,x,y,alpha,beta)

% HMTX_MVM performs the H-matrix full vector product y = alpha*H*x+beta*y
%
% USE:
% y = hmtx_mvm(H,x,y,alpha,beta)
%
% INPUTS:
% 'H': H-matrix structure
% 'x': input full vector
% 'y': input/output full vector
% 'alpha': scalar coefficient
% 'beta': scalar coefficient
%
% OUTPUTS:
% 'y': full vector such that y = alpha*H*x+beta*y
%
% NOTE:
%
% VERSION:
% Date: 21.12.2012
% Copyright(C) 2012-2014: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 08.01.2013: grouped partition of vector
% 10.01.2013: fixed bug of indexing
% 17.01.2013: generalized for matrices
% 21.01.2013: changed name from HMTX_AX
% 21.01.2013: changed oreder of operations in the case of 'rkmatrix'
% 13.02.2014: new routine with input vector
% 14.02.2014: added check on existence of matrix content
% 24.02.2014: x,y can can be full matrices

if strcmpi(H.type,'supermatrix')
    y = hmtx_mvm(H.M{1,1},x,y,alpha,beta);
    y = hmtx_mvm(H.M{2,1},x,y,alpha,beta);
    y = hmtx_mvm(H.M{1,2},x,y,alpha,beta);
    y = hmtx_mvm(H.M{2,2},x,y,alpha,beta);
elseif strcmpi(H.type,'fullmatrix') && ~isempty(H.M)
    y(H.irow,:) = alpha*H.M*x(H.jcol,:)+beta*y(H.irow,:);
elseif strcmpi(H.type,'rkmatrix') && ~isempty(H.U)
    y(H.irow,:) = alpha*H.U*(x(H.jcol,:)'*H.V)'+beta*y(H.irow,:);
end

end