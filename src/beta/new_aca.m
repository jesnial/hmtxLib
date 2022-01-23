function [U,V,k,err] = new_aca(fkern,irow,jcol,tol)
% ACA method according to Rjasanow p. 126
% computes the adaptive cross-approximation of a matrix A 
%
%INPUTS : 
% column resp. row indices ('jcol', resp. 'irow')
% 'fkern' kernel function
% 'tol' tolerance respectively to the Frobenius norm
%
%OUTPUTS :
% 'U' matrix with rows
% 'V' matrix with columns
% 'k' matrix rank (counter)
% 'err' final tolerance on Frobenius matrix norm

%
%changed from appending elements to a list to maring them with 0 or 1 for
%pivot indices
M = length(irow);
N = length(jcol);


%Initialisation
S0=0;
I=zeros(N,1);
J=zeros(1,M);
k=0;
c=zeros(N,1);
r=zeros(1,M);
kmax=min(N,M);
U = zeros(M,kmax); % left low rank matrix
V = zeros(N,kmax); % right low rank matrix
ipvt=1;
CrossType='Column';
err=1;

while err>tol && k<min(N,M)
    k=k+1
    %start with the next not yet generated row
    %     idxi=ismember(irow,I);
    %     ipvt_list=find(~idxi);
    %     ipvt=ipvt_list(1,1);
    %CrossType='Row';

    %Type of the Cross
    %if strcmp(CrossType, 'Row')
    %Generate row
    a= fkern(irow(ipvt),jcol)

    %update control vector
    I(ipvt)=1;
    r=r+abs(a);
   %end

    %test
    if abs(a)==0
        break;
    end

    %row of the residual
    rv=a-U(ipvt,1:k)*V(:,1:k)';

    %pivot column
    %idxj=ismember(jcol,J)
    [~,jpvt] = max(abs(rv.*(~J)));
    J(jpvt)=1;

    %test
    if abs(rv)~=0
        %normalizing constant
        gamma=1./rv(jpvt);

        %generate column, update control vector
        b=fkern(irow,jcol(jpvt));
        %append(J,jpvt);
        c=c+abs(b);

        %test
        if abs(b)==0
            break
        end

        %column of the residual and the pivot row
        ru=b-U(:,1:k)*V(jpvt,1:k)'
        %idxi=ismember(irow,I);
        [~,ipvt]=max(abs(ru.*(~I)));

        %test
        if abs(ru)~=0
            gamma=1./ru(ipvt)


            %generate row, uodate control vector
            a= fkern(irow(ipvt),jcol);
            %append(I,ipvt);
            I(ipvt)=1;
            r=r+abs(a);
            %row of the residual
            rv=a-U(ipvt,1:k)*V(:,1:k)'
            %pivot column
            %idxj=ismember(jcol,J)
            %[~,jpvt] = max(abs(rv.*(~idxj))); %check this, might not work... twice the same value ?
            %append(J,jpvt);

            %new vectors
            u = gamma*ru
            v = rv(:)
            % new approximation
            U(:,k) = u;
            V(:,k) = v;

            %frobenius norm of the approximation

            napprox=norm(U(:,1:k)*V(:,1:k)', 'fro')
            normunormv=norm(u,'fro')*norm(v,'fro')

            %test
            err=normunormv/napprox

            end
        end

        %check control vectors
        %ok because we go through everything in the while loop ?

        %then crop matrices

    end
  
end
        

        
       
    
    
   
    
    
    
        
        
    


