function H = hmtx_create(type,irow,jcol)

% HMTX_CREATE creates an empty H-matrix
%
% USE:
% H = hmtx_create(type,irow,jcol)
%
% INPUTS:
% 'type': H-matrix type (Available: 'supermatrix', 'rkmatrix', 'fullmatrix')
% 'irow': row indices
% 'jcol': col indices
%
% OUTPUTS:
% 'H': H-matrix structure
%
% NOTE:
%
% VERSION:
% Date: 23.01.2013
% Copyright(C) 2013-2019: Fabio Freschi (fabio.freschi@polito.it)
%                         
%
% HISTORY:
% 04.09.2014: help added
% 21.06.2018: field K added


H.type = type;            % H-matrix type
H.irow = int32(irow(:));  % global indices of rows
H.jcol = int32(jcol(:));  % global indices of columns
H.nrow = length(irow);    % number of rows
H.ncol = length(jcol);    % number of columns

% matrix data or cell array
if strcmpi(type,'supermatrix')
    H.M = cell(2,2);
else
    H.M = [];
end
H.U = [];      % ACA approximation M = U*V'
H.V = [];      % ACA approximation M = U*V'
H.eps = [];    % tolerance of ACA
H.k = [];      % actual rkmatrix rank
H.kMax = [];   % maximum rkmatrix rank


% H = struct('type',type,...
%     'irow',int32(irow(:)),...
%     'jcol',int32(jcol(:)),...
%     'nrow',length(irow),...
%     'ncol',length(jcol),...
%     'M',[],...
%     'U',[],...      % ACA approximation M = U*V'
%     'V',[],...      % ACA approximation M = U*V'
%     'eps',[],...    % tolerance of ACA
%     'k',[],...      % actual rkmatrix rank
%     'kMax',[]);     % maximum rkmatrix rank
% 
% if strcmpi(type,'supermatrix')
%     H.M = cell(2,2);
% end

% % 7.6 s

end