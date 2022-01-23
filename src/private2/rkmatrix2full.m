function Af = rkmatrix2full(Ark)

% RKMATRIX2FULL converts rkmatrix to fullmatrix
%
% USE:
% Af = rkmatrix2full(Ark)
%
% INPUTS:
% 'Ark': H-matrix in rkmatrix format
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
% 17.02.2014: compliant with fake empty rkmatrices
% 07.12.2019: separated input/output

% matrix structure
Af = hmtx_create('fullmatrix',Ark.irow,Ark.jcol);

% load fullmatrix
if Ark.k == 0
    Af.M = zeros(Ark.nrow,Ark.ncol);
else
    Af.M = Ark.U*Ark.V';
end

% tolerance
Af.eps = Ark.eps;

end
