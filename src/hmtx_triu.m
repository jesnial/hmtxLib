function U = hmtx_triu(H)

% HMTX_TRIU extracts the upper triangular part of a H-matrix
%
% USE:
% U = hmtx_triu(H)
%
% INPUTS:
% 'H': H-matrix structure
%
% OUTPUTS:
% 'L': upper triangular part of H
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
% 29.12.2019: uses sparse for empty off-diagonal blocks

if strcmpi(H.type,'supermatrix')
    U = hmtx_create('supermatrix',H.irow,H.jcol);
    
    % create fake off-diagonal rkmatrices
    U.M{2,1} = hmtx_create('rkmatrix',H.M{2,1}.irow,H.M{2,1}.jcol);
    U.M{2,1}.U = sparse([],[],[],U.M{2,1}.nrow,0);
    U.M{2,1}.V = sparse([],[],[],U.M{2,1}.ncol,0);
    U.M{2,1}.eps = 0;
    U.M{2,1}.k = 0;
    
    % copy lower block as it is
    U.M{1,2} = H.M{1,2};
    
    % block 1,1
    if strcmpi(H.M{1,1}.type,'supermatrix')
        % create diagonal blocks matrices
        U.M{1,1} = hmtx_create('supermatrix',H.M{1,1}.irow,H.M{1,1}.jcol);
        U.M{1,1} = hmtx_triu(H.M{1,1});
    elseif strcmpi(H.M{1,1}.type,'fullmatrix')
        U.M{1,1} = hmtx_create('fullmatrix',H.M{1,1}.irow,H.M{1,1}.jcol);
        U.M{1,1}.M = triu(H.M{1,1}.M);
        U.M{1,1}.eps = 0;
    else
        error('HMTX_TRIU: Diagonal block (1,1) in rkmatrix format\n');
    end
    
    % block 2,2
    if strcmpi(H.M{2,2}.type,'supermatrix')
        % create diagonal blocks matrices
        U.M{2,2} = hmtx_create('supermatrix',H.M{2,2}.irow,H.M{2,2}.jcol);
        U.M{2,2} = hmtx_triu(H.M{2,2});
    elseif strcmpi(H.M{2,2}.type,'fullmatrix')
        U.M{2,2} = hmtx_create('fullmatrix',H.M{2,2}.irow,H.M{2,2}.jcol);
        U.M{2,2}.M = triu(H.M{2,2}.M);
        U.M{2,2}.eps = 0;
    else
        error('HMTX_TRIU: Diagonal block (2,2) in rkmatrix format\n');
    end
end


end

