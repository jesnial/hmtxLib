function H = hmtx_diagonal(H,d)

% fill matrix
H = filldiagonal(H,d);

% compress 'rkmatrix' and 'fullmatrix' blocks
H = hmtx_compress(H);

end

% recursive function
function H = filldiagonal(H,d)

if strcmpi(H.type,'supermatrix')
    H.M{1,1} = filldiagonal(H.M{1,1},d);
    H.M{2,2} = filldiagonal(H.M{2,2},d);
    
    H.M{2,1} = hmtx_create('rkmatrix',H.M{2,1}.irow,H.M{2,1}.jcol);
    H.M{2,1}.U = zeros(H.M{2,1}.nrow,0);
    H.M{2,1}.V = zeros(H.M{2,1}.ncol,0);
    H.M{2,1}.k = 0;
    H.M{2,1}.eps = 0;

    H.M{1,2} = hmtx_create('rkmatrix',H.M{1,2}.irow,H.M{1,2}.jcol);
    H.M{1,2}.U = zeros(H.M{1,2}.nrow,0);
    H.M{1,2}.V = zeros(H.M{1,2}.ncol,0);
    H.M{1,2}.k = 0;
    H.M{1,2}.eps = 0;
elseif strcmpi(H.type,'fullmatrix')
    H.M = diag(d(H.irow),0);
    H.eps = 0;
elseif strcmpi(H.type,'rkmatrix')
    fprintf('oooooo\n')
end

end