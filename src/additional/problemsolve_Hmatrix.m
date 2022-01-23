function [x,xaug] = problemsolve_Hmatrix(R,Lk,b0,f,idfix,xfix,idfloat,varargin)

% PROBLEMSOLVE_HMATRIX solves the problem Ax = bsolves the problem Ax = b
% applying fixed and floating boundery conditions
%
% USE:
% x = problemsolve_Hmatrix(R,Lk,b0,f,idfix,xfix,idfloat,varargin)
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

t1 = tic;
fprintf('* system solution... ');

% default
errMAX = 1e-6;
iterMAX = 200;
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

% number of unknowns
nx = Lk.nrow;

%% floating node matrix
if ~isempty(idfloat)
    ir = [];
    jc = [];
    bx = [];
    for i = 1:length(idfloat)
        % number of floating nodes
        nfloat = length(idfloat{i});
        % filling
        ir = [ir
            (i-1)*(nfloat-1)+(1:nfloat-1).'  % master node
            (i-1)*(nfloat-1)+(1:nfloat-1).'  % slave nodes
            ];
        jc = [jc
            repmat(idfloat{i}(1),nfloat-1,1)  % master node
            idfloat{i}(2:end)                 % slave nodes
            ];
        bx = [bx
            ones(nfloat-1,1)   % master node
            -ones(nfloat-1,1)  % slave nodes
            ];
    end
    B = sparse(ir(:),jc(:),bx(:),max(ir),nx);
    % add floating boundary contribution
    nrb = size(B,1);
    
    % preallocate matrix
    Ak = hmtx_create('supermatrix',(1:Lk.nrow+nrb),(1:Lk.ncol+nrb));
    Ak.eps = Lk.eps;
    
    % block 2,1
    Ak.M{2,1} = hmtx_create('fullmatrix',(Lk.nrow+1:Lk.nrow+nrb),(1:Lk.ncol));
    Ak.M{2,1}.M = B;
    Ak.M{2,1}.eps = Lk.eps;
    % block 1,2
    Ak.M{1,2} = hmtx_create('fullmatrix',(1:Lk.nrow),(Lk.ncol+1:Lk.ncol+nrb));
    Ak.M{1,2}.M = B';
    Ak.M{1,2}.eps = Lk.eps;
    % block 2,2
    Ak.M{2,2} = hmtx_create('fullmatrix',(Lk.nrow+1:Lk.nrow+nrb),(Lk.ncol+1:Lk.ncol+nrb));
    Ak.M{2,2}.M = sparse(Ak.M{2,2}.nrow,Ak.M{2,2}.ncol);
    Ak.M{2,2}.eps = Lk.eps;
    
    %% rhs
    b0 = [b0; zeros(nrb,1)];
    x0 = [x0; zeros(nrb,1)];
else
    % preallocate matrix
    Ak = hmtx_create('supermatrix',(1:Lk.nrow),(1:Lk.ncol));
    Ak.eps = Lk.eps;
    
    % block 2,1
    Ak.M{2,1} = hmtx_create('fullmatrix',[],[]);
    Ak.M{2,1}.eps = Lk.eps;
    % block 1,2
    Ak.M{1,2} = hmtx_create('fullmatrix',[],[]);
    Ak.M{1,2}.eps = Lk.eps;
    % block 2,2
    Ak.M{2,2} = hmtx_create('fullmatrix',[],[]);
    Ak.M{2,2}.eps = Lk.eps;

end

%% resistance
[ir,jc,r] = find(R);

%% solution
% preallocation
x = zeros(nx,length(f));
for i = 1:length(f)
    % block 1,1
    Ak.M{1,1} = hmtx_tH(Lk,1i*2*pi*f(i));
    
    % add the resistance
    Ak.M{1,1} = hmtx_addentry(Ak.M{1,1},ir,jc,r);
    
    % rhs
    b = -1i*2*pi*f(i)*b0;

    % impose fixed boundary conditions
    [Ak.M{1,1},b] = hmtx_fixedboundary(Ak.M{1,1},b,idfix,xfix);
    
    % anonymous function to calculate A*x
    afun = @(x)hmtx_HxM(Ak,x);
        
    % gmres solution
    [xaug,flag,relres,iter,resvec] = gmres(afun,b,[],errMAX,iterMAX,[],[],x0);
    
    x(:,i) = xaug(1:nx);
    
    % new starting vector
    x0 = xaug;
    
    % convergence plot
    if strcmpi(gout,'y')
        figure;
        plot(1:iter(2)+1,log10(resvec/norm(b)),'r-','LineWidth',2);
        hold on
        set(gca,'FontSize',FontSize);
        plot(1:iter(2)+1,log10(errMAX)*ones(1,iter(2)+1),'b--','LineWidth',2); % convergence tolerance
        axis([0 iter(2)+1 log10(errMAX/10) 1])
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