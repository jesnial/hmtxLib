function [U,V,k,err] = aca(fKern,iRow,jCol,tol,kMax)

% ACA calculates the adaptive cross approximation of a matrix A within a
% specified tolerance such that
% A ~ U*V'
%
% USE:
% [U,V,err] = aca(fKern,iRow,jCol,tol,kMax)
%
% INPUTS:
% 'fKern': handle to kernel function  
% 'iRow','jCol': pointers to the global indexing of rows and cols of the matrix
% 'tol': tolerance for the Frobenius norm
% 'kMax': maximum rank (default: Inf: unconstrained rank)
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
% Copyright(C) 2012-2020: Fabio Freschi (fabio.freschi@polito.it)
%                         Jessie levillain (jessielevillain@gmail.com)
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
% 04.11.2019: fixed the selection of the pivoting row/col
% 04.11.2019: fixed the while condition
% 03.12.2019: fixed while condition
% 04.12.2019: fixed maximum rank calculation
% 04.12.2019: changed calculation of error
% 05.12.2019: increased robustness with exit on null rows/cols
% 03.01.2020: changed use of index k
% 10.01.2020: forced creation of rows and cols
% 10.01.2020: added kMax as input
% 13.01.2020 : fixed kMax values

% dimensions
M = length(iRow);
N = length(jCol);

% max rank allowed
kMax = min([M,N,kMax]);

% preallocation
U = zeros(M,kMax); % left low rank matrix
V = zeros(N,kMax); % right low rank matrix
I = zeros(M,1); % row pivot list
J = zeros(1,N); % column pivot list

% initialization
iPvt = 1;
I(iPvt) = 1;
err2 = Inf;
k = 0;
tol2 = tol^2;
normS2 = 0;

while err2 > tol2 && k < kMax

    % update counter
    k = k+1;        

    % generation of the row
    r = fKern(iRow(iPvt),jCol);
    % row of the residuum
    rv = r(:).'-U(iPvt,1:k)*V(:,1:k)';
    % check for null row
    if all(rv == 0)
        break;
    end
    % pivot column
    [~,jPvt] = max(abs(rv.*(~J)));
    J(jPvt) = 1;
    
    % generation of the column
    c = fKern(iRow,jCol(jPvt));
    % column of the residuum
    ru = c(:)-U(:,1:k)*V(jPvt,1:k)';
    % check for null column
    if all(ru == 0)
        break;
    end
    % pivot row
    [~,iPvt] = max(abs(ru.*(~I)));
    I(iPvt) = 1;
    
    % normalizing constant
    gamma = 1./rv(jPvt);
    % new vectors
    u = gamma*ru;
    v = rv(:);
    
    % new approximation
    U(:,k) = u;
    V(:,k) = v;

    % evaluate norm
    normu2v2 = norm(u,'fro')^2*norm(v,'fro')^2;
    % Frobenius norm
    normS2 = normS2+2*sum(real((u.'*U(:,1:k-1)).*(v.'*V(:,1:k-1))))+normu2v2;
    
    % error check
    err2 = normu2v2/normS2;

end

% remove unused rows/columns
U = U(:,1:k);
V = V(:,1:k);

% save local error 
err = sqrt(err2);

end