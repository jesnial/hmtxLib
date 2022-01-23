function c = kernelcol(jcol,P,Q)
    c = 1./distance(Q,P(jcol,:));
end