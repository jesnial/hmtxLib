function H2 = reduce_memory(H)

H2.irow = H.irow;
H2.jcol = H.jcol;

if strcmpi(H.type,'supermatrix')
    H2.M = cell(2,2);
    H2.M{1,1} = reduce_memory(H.M{1,1});
    H2.M{2,1} = reduce_memory(H.M{2,1});
    H2.M{1,2} = reduce_memory(H.M{1,2});
    H2.M{2,2} = reduce_memory(H.M{2,2});
elseif strcmpi(H.type,'fullmatrix')
    H2.M = H.M;
elseif strcmpi(H.type,'rkmatrix')
    H2.U = H.U;
    H2.V = H.V;
end
    
    
