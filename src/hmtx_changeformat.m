function Hout = hmtx_changeformat(Hin,Hout)

% HMTX_CHANGEFORMAT converts a generic H-matrix to another H-matrix with a
% different cluster tree
%
% USE:
% Hout = hmtx_changeformat(Hin,Hout)
%
% INPUTS:
% 'Hin': input H-matrix
% 'Hout': output H-matrix with different format
%
% OUTPUTS:
% 'Hout': output H-matrix filled with the netries of Hin
%
% NOTE:
%
% VERSION:
% Date: 05.01.2020
% Copyright(C) 2020: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:

if strcmpi(Hin.type,Hout.type)
    Hout = Hin;
else
    if strcmpi(Hin.type,'supermatrix') && strcmpi(Hout.type,'fullmatrix')
        Hout = super2full(Hin);
    elseif strcmpi(Hin.type,'supermatrix') && strcmpi(Hout.type,'rkmatrix')
        Hout = super2rkmatrix(Hin,Hout.kMax);
    elseif strcmpi(Hin.type,'rkmatrix') && strcmpi(Hout.type,'fullmatrix')
        Hout = rkmatrix2full(Hin);
    elseif strcmpi(Hin.type,'fullmatrix') && strcmpi(Hout.type,'rkmatrix')
        Hout = full2rkmatrix(Hin,Hout.kMax);
    elseif strcmpi(Hin.type,'rkmatrix') && strcmpi(Hout.type,'supermatrix')
        Hout = rkmatrix2super(Hin,Hout);
    elseif strcmpi(Hin.type,'fullmatrix') && strcmpi(Hout.type,'supermatrix')
        Hout = full2super(Hin,Hout);
    end
end