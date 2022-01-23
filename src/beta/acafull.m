function [U,V,k,err] = acafull(fkern,irow,jcol,eps)

% ACAFULL calculates the adaptive cross approximation of a matrix A within
% a specified tolerance such that A ~ U*V' using full pivoting
%
% USE:
% [U,V,err] = acafull(fkern,irow,jcol,eps)
%
% INPUTS:
% 'fkern': handle to kernel function  
% 'irow','jcol': pointers to the global indexing of rows and cols of the matrix
% 'eps': tolerance for the Frobenius norm
%
% OUTPUTS:
% 'U': matrix with rows
% 'V': matrix with columns
% 'k': matrix rank
% 'err': final tolerance on Frobenius matrix norm
%
% NOTE:
% Alg 3.2, S. Rjasanow, O. Steinbach, "The fast solution of Boundary Integral
% Equations", Springer, 2007
%
% VERSION:
% Date: 17.12.2012
% Copyright(C) 2012-2014: Fabio Freschi (fabio.freschi@polito.it)
% 
% HISTORY
% 09.01.2013: rewritten in function form
% 25.08.2014: cleanup

% calculate the complete matrix K
R = fkern(irow,jcol);  % here R = K
normK = norm(R,'fro');

% initialization
err = 1;
k = 0;

% dimensions
M = length(irow);
N = length(jcol);

% max rank allowed
nrmax = 20;

% preallocation
U = zeros(M,nrmax); % left low rank matrix
V = zeros(N,nrmax); % right low rank matrix

while err > eps && k <= min(M,N)
    % update counter
    k = k+1;

    % increase dimension if needed
    if k > size(U,2)
        U = [U zeros(M,nrmax)];
        V = [V zeros(N,nrmax)];
    end

    % find pivot
    [~,imaxR] = max(R(:));
    [imax,jmax] = ind2sub([M N],imaxR);
    % normalizing constant
    gamma = 1./R(imaxR);
    % new functions
    u = gamma*R(:,jmax);
    v = R(imax,:);
    % new residual
    R = R-u*v;
    % new approximation
    U(:,k) = u;
    V(:,k) = v(:);
    
    % norm
    normR = norm(R,'fro');
    err = normR/normK;
end

% remove unused rows/columns
U = U(:,1:k);
V = V(:,1:k);

end
