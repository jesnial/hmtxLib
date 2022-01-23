function h = hmtx_plot(A)

% HMTX_PLOT plots the structure of the H-matrix
%
% USE:
% h = hmtx_plot(A)
%
% INPUTS:
% 'A': H-matrix as created by HMTX_CLUSTER
%
% OUTPUTS:
% 'h': figure handle
%
% NOTE:
%
% VERSION:
% Date: 21.12.2012
% Copyright(C) 2012-2019: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 07.01.2013: exchanged x-y coordinates
% 07.01.2013: added rank
% 08.01.2013: added color maps instead of text
% 23.04.2013: plot all patches at the same time
% 06.11.2015: added DRAWNOW to force screen update
% 09.11.2019: fixed plot of a single fullmatrix
% 02.12.2019: fixed position of tick labels in colorbar

% starting row/col numbers
ir = 0;
ic = 0;
% initialize points and triangulations
Pf = [];
Tf = [];
Pk = [];
Tk = [];
Mk = [];

% call plot recursive function
[Pf,Tf,Pk,Tk,Mk] = patchhmatrix(A,ir,ic,Pf,Tf,Pk,Tk,Mk);

% create void figure
h = figure;
hold on, axis equal, grid off
axis([0 A.ncol+1 0 A.nrow+1]);

% reverse y-axis
set(gca,'YDir','reverse');

% plot patches
patch('Faces',Tf,'Vertices',Pf,'FaceColor','w','EdgeColor','k','LineWidth',1);
if size(Mk,2) == 1
    patch('Faces',Tk,'Vertices',Pk,'FaceVertexCData',Mk,'FaceColor','flat','EdgeColor','k','LineWidth',1);

    % values
    scaleValues = min(Mk):1:max(Mk); %unique(Mk);
    % colorbar
    hb = colorbar;
    colormap(parula(length(scaleValues)));
    % discrete levels
    step = (length(scaleValues)-1)/length(scaleValues);
    hb.Ticks = min(Mk)+step/2:step:max(Mk);
    hb.TickLabels = cellstr(num2str(scaleValues(:)));
    
elseif size(Mk,2) == 3
    patch('Faces',Tk,'Vertices',Pk,'FaceColor',Mk,'EdgeColor','k','LineWidth',1);
end

% setup figure

% set white background
set(h,'Color',[1 1 1]);
set(gca,'FontSize',16);

drawnow
%print('-depsc2','-r600','Htest.eps')

end

function [Pf,Tf,Pk,Tk,Mk] = patchhmatrix(A,ir,ic,Pf,Tf,Pk,Tk,Mk)

% plot frame
P = [ic+.5 ir+.5
    ic+.5 ir+A.nrow+.5
    ic+A.ncol+.5 ir+A.nrow+.5
    ic+A.ncol+.5 ir+.5
    ];
T = [1 2 3 4];

% plot square or recursively call patchmatrix
if strcmpi(A.type,'supermatrix')
    [Pf,Tf,Pk,Tk,Mk] = patchhmatrix(A.M{1,1},ir,ic,Pf,Tf,Pk,Tk,Mk);
    [Pf,Tf,Pk,Tk,Mk] = patchhmatrix(A.M{1,2},ir,ic+A.M{1,1}.ncol,Pf,Tf,Pk,Tk,Mk);
    [Pf,Tf,Pk,Tk,Mk] = patchhmatrix(A.M{2,1},ir+A.M{1,1}.nrow,ic,Pf,Tf,Pk,Tk,Mk);
    [Pf,Tf,Pk,Tk,Mk] = patchhmatrix(A.M{2,2},ir+A.M{1,1}.nrow,ic+A.M{1,1}.ncol,Pf,Tf,Pk,Tk,Mk);
elseif strcmpi(A.type,'fullmatrix')
    Tf = [Tf; T+size(Pf,1)];
    Pf = [Pf; P];
elseif strcmpi(A.type,'rkmatrix')
    Tk = [Tk; T+size(Pk,1)];
    Pk = [Pk; P];
    
    if isfield(A,'U')
        Mk = [Mk; size(A.U,2)];
    else
        Mk = [Mk; [0 1 0]];
    end
end

end