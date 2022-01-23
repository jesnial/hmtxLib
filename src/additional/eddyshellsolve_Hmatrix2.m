function [x,xaug] = eddyshellsolve_Hmatrix(R,L,b0,f,idfix,xfix,idfloat,varargin)

% EDDYSHELLSOLVE_HMATRIX solves the problem Ax = bsolves the problem Ax = b
% applying fixed and floating boundery conditions
%
% USE:
% x = eddyshellsolve_Hmatrix(R,Lk,b0,f,idfix,xfix,idfloat,varargin)
%
% INPUTS:
% 'R': resistance sparse matrix
% 'Lk': inductance H-matrix
% 'b0': I0*Lc with filamentary coils, zeros() otherwise
% 'f': vector with working frequencies
% 'idfix': vector containing indices of assigned variables
% 'xfix': vector containing assigned values
% 'idfloat': cell array with groups of equivalue unknowns
% 'varargin': optional arguments
%   * 'errMAX': tolerance of the method (default: 1e-8)
%   * 'iterMAX': maximum number of iterations (default: 500)
%   * 'x0': initial guess (default: zero vector)
%   * 'gout': flag for graphic output. Available 'y'/'n' (default: 'y')
%
% OUTPUTS
% 'x': complete vector of unknowns after solution
%
% NOTE:
%
% VERSION:
% Date: 23.04.2013
% Copyright(C) 2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 26.04.2013: changed foating boundary condition managing
% 26.04.2013: added possibility of frequency sweep
% 29.04.2013: added filamentary coil rhs
% 20.11.2013: changed name from PROBLEMSOLVE_HMATRIX
% 21.02.2014: uses new HMTX_ADD function
% 22.02.2014: uses LAGRANGIANBC to calculate the matrix for boundary conditions
% 24.02.2014: uses [D 0; 0 0] preconditioner

t1 = tic;
fprintf('* system solution... ');

% default
errMAX = 1e-6;
iterMAX = 100;
x0 = zeros(size(b0,1),1);
gout = 'y';
FontSize = 16;

% check varargin
for i = 1:2:length(varargin)
    vin = varargin{i};
    if isequal(vin,'errMAX')
        errMAX = varargin{i+1};
    elseif isequal(vin,'iterMAX')
        iterMAX = varargin{i+1};
    elseif isequal(vin,'x0')
        x0 = varargin{i+1};
    elseif isequal(vin,'gout')
        gout = varargin{i+1};
    end
end

%% floating node matrix
B = lagrangianbc(idfix,idfloat,L.nrow);

% initial guess
x0 = [x0; zeros(size(B,1),1)];

% rhs lagrangian part
b2 = [xfix; zeros(size(B,1)-length(idfix),1)];

%% resistance
[ir,jc,r] = find(R);
R = hmtx_copystruct(L);
R = hmtx_addentry(R,ir,jc,r);

%% solution
% preallocation
x = zeros(L.nrow,length(f));
for i = 1:length(f)
    % A = R+i1*w*L
    A = hmtx_add(hmtx_tH(L,1i*2*pi*f(i)),R,1,1);
    
    % impose fixed boundary conditions
    %[A,b] = hmtx_fixedboundary(A,b,idfix,xfix);
    
    % transform the system directly applying the preconditioner
    
    % factorize A
    [LL,UU] = hmtx_ssor(A);
    
    % rhs
    b = [-1i*2*pi*f(i)*hmtx_usolve(UU,hmtx_lsolve(LL,b0)); B*hmtx_usolve(UU,hmtx_lsolve(LL,B'*b2))];
    
    % anonymous function to calculate A*x
    if isempty(B)
        afun = @(x)hmtx_mvm(A,x(1:A.nrow),zeros(A.nrow,1),1,1);
    else
        afun = @(x)[hmtx_usolve(UU,hmtx_lsolve(LL,hmtx_mvm(A,x(1:A.nrow),zeros(A.nrow,1),1,1)))+hmtx_usolve(UU,hmtx_lsolve(LL,B'*x(A.nrow+1:end)));
            B*hmtx_usolve(UU,hmtx_lsolve(LL,B'*B*x(1:A.nrow)))];
    end
        
    % gmres solution
    [xaug,~,relres,~,resvec] = gmres(afun,b,[],errMAX,iterMAX,[],[],x0);
    
    x(:,i) = xaug(1:L.nrow);
    
    % new starting vector
    x0 = xaug;
    
    % convergence plot
    if strcmpi(gout,'y')
        figure;
        plot(1:length(resvec),log10(relres*resvec/resvec(end)),'r-','LineWidth',2);
        hold on
        set(gca,'FontSize',FontSize);
        plot(1:length(resvec),log10(errMAX)*ones(1,length(resvec)),'b--','LineWidth',2); % convergence tolerance
        axis([0 length(resvec) log10(errMAX/10) 1])
        xlabel('iteration #','FontSize',FontSize)
        ylabel('log10(|| Ax-b ||/|| b ||)','FontSize',FontSize);
        titstr = ['frequency: ' num2str(f(i)/1e3) ' kHz'];
        title(titstr);
        % set white background
        set(gcf,'Color',[1 1 1]);
    end
    
end

fprintf('done %6.2f sec\n',toc(t1));

end