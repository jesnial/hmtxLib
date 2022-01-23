function [x,normErr] = hmtx_backslash(A,b,varargin)

% 18.01.2019
% varargin
% * 'symmetric/hermitian'
% * 'iterative refinement'
iterMAX = 3;

% LU factors
[L,U] = hmtx_lu(A);

% solution
x = hmtx_usolve(U,hmtx_lsolve(L,b));

% preallocation
normErr = zeros(iterMAX);

% iterative refinement
for i = 1:iterMAX
    r = b - hmtx_HxM(A,x);
    e = hmtx_usolve(U,hmtx_lsolve(L,r));
    % update
    x = x+e;
    % error norm
    normErr(i) = norm(e);
end
