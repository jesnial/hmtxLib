function L = mim_Lcoeff(irow,jcol,P,F)

% 'irow': global indices of rows (points)
% 'jcol': global indices of columns (faces)
% 'P': nodes of triangulaiton
% 'F': faces of triangulation
% 'Q': field points

% 26.04.2013

% 20.11.2013: enforce triangulation to be double matrix
persistent TRI

if isempty(TRI)
  % build data structure
  TRI = TriRep(double(F),P);
end

% find triangles attached to irow
itri = vertexAttachments(TRI,double(irow(:)));
% find unique set of elements
itri = unique(cell2mat(itri.'));
% get elements
Ti = F(itri,:);
% get unique nodes
[ipnt,~,n] = unique(Ti);
% points
Pi = P(ipnt,:);
Ti = reshape(n,size(Ti));

% find triangles attached to jcol
jtri = vertexAttachments(TRI,double(jcol(:)));
% find unique set of elements
jtri = unique(cell2mat(jtri.'));
% get elements
Tj = F(jtri,:);
% get unique nodes
[jpnt,~,n] = unique(Tj);
% points
Pj = P(jpnt,:);
Tj = reshape(n,size(Tj));

try
    L = mutualinductancematrixf(Pi,int32(Ti),Pj,int32(Tj));
catch
    L = mutualinductancematrixm(Pi,Ti,Pj,Tj);
end

% find irow in local numbering
[~, loc] = ismember(ipnt,irow);
irowloc(loc(loc ~= 0)) = find(loc);
% find jcol in local numbering
[~, loc] = ismember(jpnt,jcol);
jcolloc(loc(loc ~= 0)) = find(loc);

L = L(irowloc,jcolloc);

end