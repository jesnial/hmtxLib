function TF = ishmtx(A)

% ISHMTX check if the input is a valid H-matrix 
%
% USE:
% TF = ishmtx(A)
%
% INPUTS:
% 'A': full matrix
% 'H': H-matrix structure
%
% OUTPUTS:
% 'TF': logical 1 (true) or 0 (false)
%
% NOTE:
%
% VERSION:
% Date: 04.09.2014
% Copyright(C) 2014-2020: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 21.06.2018: added K field
% 03.01.2020: added KMAX field

TF = isstruct(A) && isfield(A,'irow') && isfield(A,'jcol') && ...
    isfield(A,'nrow') && isfield(A,'ncol') && ...
    isfield(A,'M') && isfield(A,'U') && isfield(A,'V') && ...
    isfield(A,'eps') && isfield(A,'k') && isfield(A,'kMax');