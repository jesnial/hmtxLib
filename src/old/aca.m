function [U,V,k,err] = aca(fkern,irow,jcol,tol)

% ACA calculates the adaptive cross approximation of a matrix A within a
% specified tolerance such that
% A ~ U*V'
%
% USE:
% [U,V,err] = aca(fkern,irow,jcol,tol)
%
% INPUTS:
% 'fkern': handle to kernel function  
% 'irow','jcol': pointers to the global indexing of rows and cols of the matrix
% 'tol': tolerance for the Frobenius norm
%
% OUTPUTS:
% 'U': matrix with rows
% 'V': matrix with columns
% 'k': matrix rank
% 'err': final tolerance on Frobenius matrix norm
%
% NOTE:
% See S. Rjasanow, O. Steinbach "The Fast Solution of Boundary Integral
% Equations Algorithm" 3.9 pp. 126--129
%
% VERSION:
% Date: 20.12.2012
% Copyright(C) 2012-2014: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 05.01.2013: FKERN as unique handle function
% 07.01.2013: IROW and JCOL passed as inputs
% 08.01.2013: added condition when pivot is chosen randomly
% 13.01.2013: creates V and not V'
% 18.01.2013: added rank to outputs
% 01.02.2013: correction on dimension of J and rv
% 16.04.2013: U and V are reallocated when necessary
% 22.04.2013: added additional exit condition to WHILE loop
% 15.07.2014: fixed bug in norm calculation
% 22.07.2014: fixed bug in the calculation of the rank
% 04.09.2014: change in the norm calculation

% dimensions
M = length(irow);
N = length(jcol);

% max rank allowed
nrmax = 20;

% preallocation
U = zeros(M,nrmax); % left low rank matrix
V = zeros(N,nrmax); % right low rank matrix
I = zeros(M,1); % row pivot list
J = zeros(1,N); % column pivot list

% initialization
ipvt = 1;
I(ipvt) = 1;
err2 = 1e6;
k = 0;
tol2 = tol^2;

while err2 > tol2 && k <= min(M,N)
    % update counter
    k = k+1;

    % increase dimension as needed
    if k > size(U,2)
        U = [U zeros(M,nrmax)];
        V = [V zeros(N,nrmax)];
    end
    
    % generation of the row
    r = fkern(irow(ipvt),jcol);
    % row of the residuum
    rv = r-U(ipvt,1:k)*V(:,1:k)';
    % check for null row
    if all(rv == 0)
        error('rv = 0');
    end
    % pivot column
    [~,jpvt] = max(abs(rv.*(~J)));
    J(jpvt) = 1;
    
    % generation of the column
    c = fkern(irow,jcol(jpvt));
    % column of the residuum
    ru = c-U(:,1:k)*V(jpvt,1:k)';
    % check for null column
    if all(ru == 0)
        error('ru = 0');
    end
    
    % pivot row
    [~,ipvt] = max(abs(ru.*(~I)));
    % check row index
    if I(ipvt) == 1
        ipvt = find(I == 0 & ru ~= 0,1,'first');
    end
    I(ipvt) = 1;
    
    % normalizing constant
    gamma = 1./rv(jpvt);
    % new vectors
    u = gamma*ru;
    v = rv(:);
    % new approximation
    U(:,k) = u;
    V(:,k) = v;

    % evaluate norm
    normu2v2 = norm(u)^2*norm(v)^2;
    % Frobenius norm of a product as trace (square)
    normS2 = trace((U(:,1:k)'*U(:,1:k))*(V(:,1:k)'*V(:,1:k)));
    % error check
    err2 = normu2v2/real(normS2);

%     % recusive norm
%     normu2v2 = norm(u)^2*norm(v)^2;
%     uuvv = 0;
%     for j = 1:k-1
%         uuvv = uuvv+2*(U(:,j)'*U(:,k))*(V(:,j)'*V(:,k));
%     end
%     normS2 = normS2+uuvv+normu2v2;
% 
%     err2 = normu2v2/normS2;
    
        
end

% remove unused rows/columns
U = U(:,1:k);
V = V(:,1:k);

% save local error 
err = sqrt(err2);

end