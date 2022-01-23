function [coeff] = fun_my_kernel(ii,jj,P,eeps)
Nii=length(ii);
Njj=length(jj);
coeff=zeros(Nii,Njj);
for a = 1:Njj
    coeff(:,a)=1./((fun_my_norm(P(ii,:)-P(jj(a),:))+eeps));
end
end

