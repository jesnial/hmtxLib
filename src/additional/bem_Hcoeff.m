function H = bem_Hcoeff(irow,jcol,P,F,Q)

% 'irow': globel indices of rows (points)
% 'jcol': globel indices of columns (faces)
% 'P': nodes of triangulaiton
% 'F': faces of triangulation
% 'Q': field points

try
    [~,H] = bem_matrix3df(P,int32(F(jcol,:)),Q(irow,:));
catch
    [~,H] = bem_matrix3dm(P,F(jcol,:),Q(irow,:));
end
% % replace diagonal elements with 1/2
% [~,rowSub,colSub] = intersect(irow,jcol);
% H(sub2ind(size(H),rowSub,colSub)) = .5;

end