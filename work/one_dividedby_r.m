function D = one_dividedby_r(irow,jcol,P,Q)

alpha = 1;

D = 1./(alpha+distance(Q(irow,:),P(jcol,:)));

end