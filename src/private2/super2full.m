function Af = super2full(As)

% SUPER2FULL converts a supermatrix into a fullmatrix
%
% USE:
% Af = super2full(As)
%
% INPUTS:
% 'As': H-matrix in supermatrix format
%
% OUTPUTS:
% 'Af': H-matrix in fullmatrix format
%
% NOTE:
%
% VERSION:
% Date: 06.02.2014
% Copyright(C) 2014-2019: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 07.12.2019: separated input/output
% 13.12.2019: added eps field

if strcmpi(As.type,'supermatrix')
    Af = hmtx_create('supermatrix',As.irow,As.jcol);
    % convert each block into rkmatrix
    Af.M{1,1} = super2full(As.M{1,1});
    Af.M{2,1} = super2full(As.M{2,1});
    Af.M{1,2} = super2full(As.M{1,2});
    Af.M{2,2} = super2full(As.M{2,2});
    % compress the 4 fullmatrix blocks
    Af = hmtx_compress(Af);
elseif strcmpi(As.type,'rkmatrix')
    Af = rkmatrix2full(As);
elseif strcmpi(As.type,'fullmatrix')
    Af = As;
end
Af.eps = As.eps;

end