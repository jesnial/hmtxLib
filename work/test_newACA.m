clear variables, close all

P = [0 0 0; 1 1 0; 1 0 0; 1 -1 0];
T = [1 3 2; 1 4 3];
M = 2*ones(size(T,1),1);

primal = createprimal2d(P,T,M);
mat = setmatproperties(M);

L = inductancematrix(primal,mat)
inod = [1 2 3 4];
clear mim_Lcoeff
L2 = mim_Lcoeff(inod,inod,primal)
