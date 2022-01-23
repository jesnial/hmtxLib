function Ark = super2rkmatrix(As,kMax)

% SUPER2RKMATRIX hiearchical approximation of a supermatrix A to a Rk
% matrix
%
% USE:
% Ark = super2rkmatrix(As,kMax)
%
% INPUTS:
% 'As': H-matrix in supermatrix format
% 'kMax': maximum rank
%
% OUTPUTS:
% 'Ark': H-matrix in rkmatrix format
%
% NOTE:
% See S. Borm, L. Grasedyck, W. Hackbusch, "Hierarchical Matrices", 2003,
% (revised version: June 2006), pp. 114
%
% VERSION:
% Date: 04.02.2014
% Copyright(C) 2014-2020: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 07.12.2019: separated input/output
% 13.12.2019: added eps field
% 03.01.2020: added KMAX as input

if strcmpi(As.type,'supermatrix')
    Ark = hmtx_create('supermatrix',As.irow,As.jcol);
    % convert each block into rkmatrix
    Ark.M{1,1} = super2rkmatrix(As.M{1,1},kMax);
    Ark.M{2,1} = super2rkmatrix(As.M{2,1},kMax);
    Ark.M{1,2} = super2rkmatrix(As.M{1,2},kMax);
    Ark.M{2,2} = super2rkmatrix(As.M{2,2},kMax);
    Ark.eps = As.eps;
    Ark.kMax = kMax;
    % compress the 4 rkmatrix blocks
    Ark = hmtx_compress(Ark);
elseif strcmpi(As.type,'fullmatrix')
    Ark = full2rkmatrix(As,kMax);
elseif strcmpi(As.type,'rkmatrix')
    Ark = As;
end
Ark.eps = As.eps;

end