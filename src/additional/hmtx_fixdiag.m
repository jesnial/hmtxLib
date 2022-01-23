function A = hmtx_fixdiag(A)

% HMTX_FIXDIAG replaces the diagonal of an H-matrix with 0.5


% first diagonal block
if strcmpi(A.M{1,1}.type,'supermatrix')
    A.M{1,1} = hmtx_fixdiag(A.M{1,1});
elseif strcmpi(A.M{1,1}.type,'fullmatrix')
    A.M{1,1}.M(1:(A.M{1,1}.nrow+1):A.M{1,1}.nrow^2) = .5;
end

% second diagonal block
if strcmpi(A.M{2,2}.type,'supermatrix')
    A.M{2,2} = hmtx_fixdiag(A.M{2,2});
elseif strcmpi(A.M{2,2}.type,'fullmatrix')
    A.M{2,2}.M(1:(A.M{2,2}.nrow+1):A.M{2,2}.nrow^2) = .5;
end

end