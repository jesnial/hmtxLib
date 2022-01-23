function MB = hmtx_Mbyte(H)

% HMTX_Mbyte returns the megabytes used by an H-matrix
%
% USE:
% MB = hmtx_Mbyte(H)
%
% INPUTS:
% 'H': H-matrix
%
% OUTPUTS:
% 'MB': memory usage in MB
%
% NOTE:
%
% VERSION:
% Date: 16.04.2013
% Copyright(C) 2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 17.04.2013: uses GETFIELD

MB = getfield(whos('H'),'bytes')/2^20;

end