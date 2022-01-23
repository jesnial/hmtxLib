clear variables, close all
a = (1:2000000)';
b = (1000000:3000000)';

tic
logUA = ~ismembc(a,b);
c = a(logUA);
ia = find(logUA);
toc

tic
[cc,iiaa] = setdiff(a,b);
toc
