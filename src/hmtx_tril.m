function L = hmtx_tril(H)

% HMTX_TRIL extracts the lower triangular part of a H-matrix
%
% USE:
% L = hmtx_tril(H)
%
% INPUTS:
% 'H': H-matrix structure
%
% OUTPUTS:
% 'L': lower triangular part of H
%
% NOTE:
%
% VERSION:
% Date: 16.02.2014
% Copyright(C) 2014: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 17.02.2014: some cleaning of the routine
% 22.11.2019: fixed empty rkmatrix block
% 22.11.2019: fixed fullmatrix block
% 16.12.2019: added tril of 1-level fullmatrix
% 29.12.2019: uses sparse for empty off-diagonal blocks

if strcmpi(H.type,'supermatrix')
    L = hmtx_create('supermatrix',H.irow,H.jcol);
    
    % create empty off-diagonal rkmatrix
    L.M{1,2} = hmtx_create('rkmatrix',H.M{1,2}.irow,H.M{1,2}.jcol);
    L.M{1,2}.U = sparse([],[],[],L.M{1,2}.nrow,0);
    L.M{1,2}.V = sparse([],[],[],L.M{1,2}.ncol,0);
    L.M{1,2}.eps = 0;
    L.M{1,2}.k = 0;
    
    % copy lower block as it is
    L.M{2,1} = H.M{2,1};
    
    % block 1,1
    if strcmpi(H.M{1,1}.type,'supermatrix')
        % create diagonal blocks matrices
        L.M{1,1} = hmtx_create('supermatrix',H.M{1,1}.irow,H.M{1,1}.jcol);
        L.M{1,1} = hmtx_tril(H.M{1,1});
    elseif strcmpi(H.M{1,1}.type,'fullmatrix')
        L.M{1,1} = hmtx_create('fullmatrix',H.M{1,1}.irow,H.M{1,1}.jcol);
        L.M{1,1}.M = tril(H.M{1,1}.M);
        L.M{1,1}.eps = H.M{1,1}.eps;
    else
        error('HMTX_TRIL: Diagonal block (1,1) in rkmatrix format\n');
    end
    
    % block 2,2
    if strcmpi(H.M{2,2}.type,'supermatrix')
        % create diagonal blocks matrices
        L.M{2,2} = hmtx_create('supermatrix',H.M{2,2}.irow,H.M{2,2}.jcol);
        L.M{2,2} = hmtx_tril(H.M{2,2});
    elseif strcmpi(H.M{2,2}.type,'fullmatrix')
        L.M{2,2} = hmtx_create('fullmatrix',H.M{2,2}.irow,H.M{2,2}.jcol);
        L.M{2,2}.M = tril(H.M{2,2}.M);
        L.M{2,2}.eps = H.M{2,2}.eps;
    else
        error('HMTX_TRIL: Diagonal block (2,2) in rkmatrix format\n');
    end
elseif strcmpi(H.type,'fullmatrix')
    L = hmtx_create('fullmatrix',H.irow,H.jcol);
    L.M = tril(H.M);
    L.eps = H.eps;
end

end

