function H = hmtx_tH(t,H)

% HMTX_TH multiplies a H-matrix by a scalar t
%
% USE:
% H = hmtx_tH(H,t)
%
% INPUTS:
% 'H': H-matrix, as created by HMTX_CLUSTER and filled by HMTX_FILL
% 't': scalar value
%
% OUTPUTS:
% 'H': H-matrix such that H = t*H
%
% NOTE:
%
% VERSION:
% Date: 22.01.2013
% Copyright(C) 2013-2020: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 03.01.2020: added help

% plot square or recursively call patchmatrix
if strcmpi(H.type,'supermatrix')
    H.M{1,1} = hmtx_tH(t,H.M{1,1});
    H.M{2,1} = hmtx_tH(t,H.M{2,1});
    H.M{1,2} = hmtx_tH(t,H.M{1,2});
    H.M{2,2} = hmtx_tH(t,H.M{2,2});    
elseif strcmpi(H.type,'fullmatrix')
    H.M = t*H.M;
elseif strcmpi(H.type,'rkmatrix')
    H.U = t*H.U;
end

end