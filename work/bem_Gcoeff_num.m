function G = bem_Gcoeff_num(irow,jcol,P,F,Q)

% 'irow': globel indices of rows (points)
% 'jcol': globel indices of columns (faces)
% 'P': nodes of triangulaiton
% 'F': faces of triangulation
% 'Q': field points

try
    [G,~] = bem_matrix3d_(P.',int32(F(jcol,:)).',Q(irow,:).');
catch
    G = bem_numeric_matrix3dm(P,F(jcol,:),Q(irow,:));
end


end