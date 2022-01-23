function h = plothmatrix(A)

% 21.12.2012

% 07.01.2013: exchanged x-y coordinates
% 07.01.2012: added rank
% 08.01.2012: added color maps instead of text

h = figure;
hold on
axis equal
axis([0 A.n+1 0 A.m+1])
set(gca,'YDir','reverse');
ir = 0;
ic = 0;
h = patchhmatrix(A,ir,ic,h);
figuresetup(h);
grid off

end

function h = patchhmatrix(A,ir,ic,h)

% plot frame
P = [ic+.5 ir+.5 
    ic+.5 ir+A.m+.5 
    ic+A.n+.5 ir+A.m+.5 
    ic+A.n+.5 ir+.5 
    ];
F = [1 2 3 4];
figure(h);
if strcmpi(A.type,'supermatrix')
    % patch('Faces',F,'Vertices',P,'FaceColor','none','EdgeColor','k');
    % plot 2x2 blocks
    h = patchhmatrix(A.M{1,1},ir,ic,h);
    h = patchhmatrix(A.M{1,2},ir,ic+A.M{1,1}.n,h);
    h = patchhmatrix(A.M{2,1},ir+A.M{1,1}.m,ic,h);
    h = patchhmatrix(A.M{2,2},ir+A.M{1,1}.m,ic+A.M{1,1}.n,h);
elseif strcmpi(A.type,'fullmatrix')
    patch('Faces',F,'Vertices',P,'FaceColor','k','EdgeColor','k','LineWidth',1);
elseif strcmpi(A.type,'rkmatrix')
    patch('Faces',F,'Vertices',P,'FaceVertexCData',size(A.U,2),'FaceColor','flat','EdgeColor','k','LineWidth',1);
end

end