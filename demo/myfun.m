function B = myfun(B)
if B.V > B.k
    B = myfun(B);
else
    B.U = 2*B.U;
    fprintf('ale\n')
end
end