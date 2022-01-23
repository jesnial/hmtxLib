clear 
close all
clc
%% carico la libreria
% dad=pwd;
% cd hm-toolbox-master\
% addpath(pwd);
% cd(dad)
%% num incognite 
Ndofs=1000;
%% creo i punti
P=zeros(Ndofs,3);
P(:,1)=linspace(0,100,Ndofs).';
%% fun kernel
eeps=1e-2;
fun=@(ii,jj)  fun_my_kernel(ii,jj,P,eeps);
%% creo sol
rng(1)
sol=rand(Ndofs,1);
%% creo rhs esatta al volo
rhs=zeros(Ndofs,1);
for ii = 1:Ndofs
    rhs(ii)=fun(ii,1:Ndofs)*sol;
end
%% opzioni di compressione per hodlr
hodlroption('block-size',20);
hodlroption('compression','svd'); % se usi function handle in realt� usa aca
hodlroption('threshold',1e-5);
%% opzioni di compressione per hss
hssoption('block-size',20);
hssoption('compression','svd'); % se usi function handle in realt� usa aca
hssoption('threshold',1e-5);
%% opzioni di compressione per hmatrix
hmatrixoption('block-size',20);
hmatrixoption('compression','svd'); % se usi function handle in realt� usa aca
hmatrixoption('threshold',1e-5);
%% HODLR
disp('====================================================================')
disp('making hodlr...')
tic
Ahodlr=hodlr('handle',fun,Ndofs,Ndofs);
toc
compr=100*getSize(Ahodlr)/(Ndofs*getSize(fun(1,1:Ndofs)));
disp(['Compression ratio hodlr=',num2str(compr),'%'])
figure 
spy(Ahodlr)
title('holdr (si pu� anche imporre una diversa divisione)')
%% HSS
disp('====================================================================')
disp('making hss from hodlr...')
Ahss=hodlr2hss(Ahodlr); % ci sono molti metodi
compr=100*getSize(Ahss)/(Ndofs*getSize(fun(1,1:Ndofs)));
disp(['Compression ratio hss=',num2str(compr),'%'])
figure
spy(Ahss)
title('hss')
%% HMATRIX
disp('====================================================================')
disp('making hmatrix...')
tic
Ahmatrix=hmatrix('handle',fun,Ndofs,Ndofs); % chiamandola cos� usa un albero di default
toc
compr=100*getSize(Ahmatrix)/(Ndofs*getSize(fun(1,1:Ndofs)));
disp(['Compression ratio hmatrix=',num2str(compr),'%'])
figure
spy(Ahmatrix)
title('hmatrix')
%% RISOLVO HODLR (\=lu)
disp('====================================================================')
disp('solving hodlr with \...')
tic
xhodlr=Ahodlr\rhs;
toc
err=100*norm(xhodlr-sol)/norm(sol);
disp(['Solution perc. error hodlr =',num2str(err),'%'])
%% RISOLVO HSS (\=lu)
disp('====================================================================')
disp('solving hss with \...')
tic
xhss=Ahss\rhs;
toc
err=100*norm(xhss-sol)/norm(sol);
disp(['Solution perc. error hss =',num2str(err),'%'])
%% RISOLVO HMATRIX (\=lu)
disp('====================================================================')
disp('solving hmatrix with \...')
tic
xhmatrix=Ahmatrix\rhs;
toc
err=100*norm(xhmatrix-sol)/norm(sol);
disp(['Solution perc. error hmatrix =',num2str(err),'%'])
%% creating cluster tree
disp('====================================================================')
disp('creating user defined cluster tree for hmatrix  ...')
sizeLeaf=30;
leaf(1,:)=0:sizeLeaf:Ndofs-sizeLeaf; %create regularly spaced vector using sizeLeaf as an increment
for ii = 1:length(leaf)-1
    leaf(2,ii)=leaf(1,ii+1)-leaf(1,ii);
end
leaf(2,end)=Ndofs-leaf(1,ii+1);
N_leaf=size(leaf,2);
eta=0.8;
%[HT] = fun_my_Hmatrix_T(N_leaf,leaf,P.',eta);
HT = hmatrix_cluster(P);
%% HMATRIX with user cluster tree
disp('====================================================================')
disp('making hmatrix with user defined cluster tree...')
hmatrixoption('threshold',1e-4); % puoi settare la tol dove vuoi
tic
Ahmatrix2=hmatrix( 'cluster', HT,'handle',fun,Ndofs,Ndofs); %
toc
compr=100*getSize(Ahmatrix2)/(Ndofs*getSize(fun(1,1:Ndofs)));
disp(['Compression ratio hmatrix2=',num2str(compr),'%'])
figure
spy(Ahmatrix2)
title('hmatrix with user cluster-tree (geometrically based)')
%% solving with gmres
disp('====================================================================')
disp('solving (lu+gmres)...')
disp('lu...')
hmatrixoption('threshold',1e-1); % cambio tol per la lu
tic
[l,u]=lu(Ahmatrix2);
toc
compr=100*(getSize(l)+getSize(u))/(Ndofs*getSize(fun(1,1:Ndofs)));
disp(['memory ratio l+u/A =',num2str(compr),'%'])
disp('gmres...')
tic
[xhmatrix2,flag,relres,iter,resvec]=my_gmres(Ahmatrix2,rhs,Ndofs,1e-12,100,l,u);
toc
err=100*norm(xhmatrix2-sol)/norm(sol);
disp(['Solution perc. error hmatrix2 gmres =',num2str(err),'%'])

