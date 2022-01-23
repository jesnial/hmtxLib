clear variables, close all


%% parameters
r0 = 1.5;    % sphere radius
h = .025;      % edge length

% surface mesh of a sphere
[P,F] = triSphere(h);

% normalize radius to r0
r = sqrt(sum(P.^2,2));
P = r0*P./repmat(r,1,3);

fid = fopen(strcat('sphere',num2str(size(F,1)),'.txt'),'w');
fprintf(fid,'%d\n',size(P,1));
fprintf(fid,'%e %e %e\n',P.');
fprintf(fid,'%d\n',size(F,1));
fprintf(fid,'%d %d %d\n',F.');
fclose(fid);

M = 2*ones(size(F,1),1);

primal = createprimal2d(P,F,M);
plotscalar2d(primal,'mesh','y');



