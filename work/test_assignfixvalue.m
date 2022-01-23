clear variables, close all

% matrix generation
alpha = 1e-2;
N = 24;
tol = 1e-8;  % tolerance of the ACA method
eta = 3;
Nmin = 10;

% create a grid of N elements in [-1,1]. The clustering algorithm works
% with 3d points
P = [linspace(-1,1,N).' zeros(N,2)];

% create the matrix 1/r, protecting the diagonal with alpha
K = 1./(alpha+distance(P,P));

% the kernel function is a simple acces to the matrix entries
fkern = @(irow,jcol)K(irow,jcol);

% create the H-matrix structure
H = hmtx_cluster(P,'eta',eta,'Nmin',Nmin);

% plot the empty H-matrix
hmtx_plot(H);

% populate H
H = hmtx_fill(H,tol,fkern);

% plot the H-matrix with the reduced ranks
hmtx_plot(H);

% rhs
b = zeros(N,1);

% idfix 
Nfix = ceil(N/10);
idfix = sort(randi([0 N],Nfix,1));
xfix = zeros(Nfix,1);

% create fullmatrix from H-matrix (only for small matrices)
Hf = getfield(hmtx_full(H),'M');

%%%%%%%%%%%%%%%%%%%%
% test: multiply by identity matrix
% create identity matrix
I = hmtx_copystruct(H);
I = hmtx_addentry(I,I.irow,I.jcol,ones(I.nrow,1));

% delete cols
H = hmtx_mult(H,I);

% Hr = hmtx_assignfixvalue(H,b,idfix,xfix);

