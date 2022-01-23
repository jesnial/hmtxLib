% ACA with partial pivoting: reconstruction of potential

clear variables, close all

% number of sources
N = 1000;
% number of targets
M = 1000;

% radius of sources
R0 = 1;
% inner radius of targets
R1 = 0;
% outer radius of targets
R2 = 10;

% desired tolerance
eps = logspace(-1,-10,10);

% source points
r = R0*rand(N,1);
theta = pi*rand(N,1);
phi = 2*pi*rand(N,1);
P = [r.*sin(theta).*cos(phi) r.*sin(theta).*sin(phi) r.*cos(theta)];

% charges
q = ones(N,1);

% target points
s = sqrt(3)*linspace(R1,max(R2),M)';
Q = [s s s];

% plot
h1 = figure; hold on, axis equal, view([1 1 1])
plot3(P(:,1),P(:,2),P(:,3),'o','MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor','b');
plot3(Q(:,1),Q(:,2),Q(:,3),'o','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r');

% calculate 1/r matrix (assume point charge q = 1)
D = distance(Q,P);
% kernel matrix K
K = 1./D;

% exact potential
tic
Vex = K*q;
toc
% norm of K
normK = norm(K,'fro')^2;
%plot(1:M,Vex,'r')

nrmax = max(round(min([M,N])/10),100);
for i = 1:length(eps)
    % initialization
    err = 1;
    imax = 1;
    
    U = zeros(M,nrmax); % left low rank matrix
    V = zeros(nrmax,N); % right low rank matrix
    I = zeros(M,1); % row pivot list
    I(imax) = 1;
    J = zeros(N,1); % column pivot list
    normS2 = 0;
    
    for k = 1:nrmax
        % generation of the row
        a = 1./distance(Q(imax,:),P);
        % row of the residuum
        rv = a-U(imax,1:k)*V(1:k,:);
        % check for null row
        if all(rv == 0)
            error('rv = 0');
        end
        % pivot column
        [~,jmax] = max(abs(rv));
        % check column index
        if J(jmax) == 1
            jmax = find(J == 0,1,'first');
        end
        J(jmax) = 1;

        % normalizing constant
        gamma = 1./rv(jmax);
        % generation of the column
        a = 1./distance(Q,P(jmax,:));
        % column of the residuum
        ru = a-U(:,1:k)*V(1:k,jmax);
        if all(ru == 0)
            error('ru = 0');
        end

        % pivot row
        [Rmax,imax] = max(abs(ru));
        % check row index
        if I(imax) == 1
            imax = find(I == 0,1,'first');
        end
        I(imax) = 1;

        % new vectors
        u = gamma*ru;
        v = rv;
        % new approximation
        U(:,k) = u;
        V(k,:) = v;
        
        % exit criterion
        normu2v2 = norm(u)^2*norm(v)^2;
        
        uuvv = 0;
        for j = 1:k-1
            uuvv = uuvv+(U(:,j)'*U(:,k))*(V(j,:)*V(k,:)');
        end
        normS2 = normS2+2*uuvv+normu2v2;
        if normu2v2/normS2 <= eps(i)^2
            break
        end
    end
    
    nr(i) = k;
    tic
    Vr(:,i) = U(:,1:k)*(V(1:k,:)*q);
    toc
end
h2 = figure; hold on
plot(sqrt(sum(Q.^2,2)),log10(abs(Vr-repmat(Vex,1,length(eps)))./repmat(Vex,1,length(eps))))
xlabel('curvilinear coordinate');
ylabel('log10 of relative error on potential');
legendCell = cellstr(num2str([eps' nr'], 'eps = %5.0e, r = %-d'));
legend(legendCell);
figuresetup(h2);
