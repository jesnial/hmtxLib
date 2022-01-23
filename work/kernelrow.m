function r = kernelrow(irow,P,Q)
    r = 1./distance(Q(irow,:),P);
end