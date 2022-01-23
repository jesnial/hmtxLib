function G = bem_Gcoeff(irow,jcol,P,F,Q)

% 'irow': global indices of rows (points)
% 'jcol': global indices of columns (faces)
% 'P': nodes of triangulaiton
% 'F': faces of triangulation
% 'Q': field points

try
    [G,~] = bem_matrix3df(P,int32(F(jcol,:)),Q(irow,:));
catch
    [G,~] = bem_matrix3dm(P,F(jcol,:),Q(irow,:));
end


end