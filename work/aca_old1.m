function [U,VT] = aca(fkern,M,N,eps)

% USE:
% [U,VT] = aca(fkern,M,N,eps)
%
% INPUTS:
% 'fkern': handle to kernel function  
% 'M','N': size of the matrix
% 'eps': tolerance for the Frobenius norm
%
% OUTPUTS:
%
% NOTE:
%
% VERSION:
% Date: 20.12.2012
% Copyright(C) 2012: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 05.01.2013: FKERN as unique handle function

% max rank allowed
nrmax = round(min([M,N])/10);

% preallocation
U = zeros(M,nrmax); % left low rank matrix
VT = zeros(nrmax,N); % right low rank matrix
I = zeros(M,1); % row pivot list
J = zeros(N,1); % column pivot list

% initialization
ipvt = 1;
I(ipvt) = 1;
normS2 = 0;

for k = 1:nrmax
    % generation of the row
    r = kernelrow(fkern,ipvt,N);
    % row of the residuum
    rv = r-U(ipvt,1:k)*VT(1:k,:);
    % check for null row
    if all(rv == 0)
        error('rv = 0');
    end
    % pivot column
    [~,jpvt] = max(abs(rv));
    % check column index
    if J(jpvt) == 1
        jpvt = find(J == 0,1,'first');
    end
    J(jpvt) = 1;
    
    % normalizing constant
    gamma = 1./rv(jpvt);
    % generation of the column
    c = kernelcol(fkern,jpvt,M);
    % column of the residuum
    ru = c-U(:,1:k)*VT(1:k,jpvt);
    if all(ru == 0)
        error('ru = 0');
    end
    
    % pivot row
    [~,ipvt] = max(abs(ru));
    % check row index
    if I(ipvt) == 1
        ipvt = find(I == 0,1,'first');
    end
    I(ipvt) = 1;
    
    % new vectors
    u = gamma*ru;
    v = rv;
    % new approximation
    U(:,k) = u;
    VT(k,:) = v;
        
    % recusive norm
    normu2v2 = norm(u)^2*norm(v)^2;
    uuvv = 0;
    for j = 1:k-1
        uuvv = uuvv+(U(:,j)'*U(:,k))*(VT(j,:)*VT(k,:)');
    end
    normS2 = normS2+2*uuvv+normu2v2;

    % exit criterion
    if normu2v2/normS2 <= eps^2
        break
    end
end

% remove unused rows/columns
U = U(:,1:k);
VT = VT(1:k,:);

end


function r = kernelrow(fkern,irow,N)

jcol = 1:N;
r = fkern(irow,jcol);

end

function c = kernelcol(fkern,jcol,M)

irow = 1:M;
c = fkern(irow,jcol);

end
