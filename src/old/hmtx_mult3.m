function C = hmtx_mult2(A,B,C)
    
% C = A*B

% 21.01.2013

% 23.01.2013: uses HMTX_CREATE to create H-matrix structure

% check dimenisons
%
% HISTORY :
%
% 29.11.2019 : fixed cross-type cases. supermatrix and rk becomes rk
% 03.12.2019 : added recompression

% convert A to C format
Ac = hmtx_copystruct(C);
Ac.eps = A.eps;
Ac = formatconversion(A,Ac);

% convert B to C format
Bc = hmtx_copystruct(C);
Bc.eps = B.eps;
Bc = formatconversion(B,Bc);

if strcmpi(C.type,'supermatrix')  
    for i = 1:2
        for j = 1:2
            fprintf('%d %d\n',i,j);
            % first term
            Ca = hmtx_copystruct(C.M{i,j});
            Ca = hmtx_mult3(Ac.M{i,1},Bc.M{1,j},Ca);

            % second term
            Cb = hmtx_copystruct(C.M{i,j});
            Cb = hmtx_mult3(Ac.M{i,2},Bc.M{2,j},Cb);

%             % first term
%             Ca = hmtx_copystruct(C.M{i,j});
% 
%             % convert to target format
%             AcMij = hmtx_copystruct(C.M{i,j});
%             AcMij.eps = Ac.M{i,1}.eps;
%             AcMij = formatconversion(Ac.M{i,1},AcMij);
%             % convert to target format
%             BcMij = hmtx_copystruct(C.M{i,j});
%             BcMij.eps = Bc.M{i,1}.eps;
%             BcMij = formatconversion(Bc.M{1,j},BcMij);
%             % multiply
%             Ca = hmtx_mult3(AcMij,BcMij,Ca);
% 
%             % second term
%             Cb = hmtx_copystruct(C.M{i,j});
% 
%             % convert to target format
%             AcMij = hmtx_copystruct(C.M{i,j});
%             AcMij.eps = Ac.M{i,2}.eps;
%             AcMij = formatconversion(Ac.M{i,2},AcMij);
%             % convert to target format
%             BcMij = hmtx_copystruct(C.M{i,j});
%             BcMij.eps = Bc.M{i,2}.eps;
%             BcMij = formatconversion(Bc.M{2,j},BcMij);
%             % multiply
%             Cb = hmtx_mult3(AcMij,BcMij,Cb);

            % sum and tolerance
            C.M{i,j} = hmtx_add(Ca,Cb,'+');
            C.M{i,j}.eps = max([Ac.M{i,1}.eps,Bc.M{1,j}.eps,Ac.M{i,2}.eps,Bc.M{2,j}.eps]);
        end
    end
elseif strcmpi(C.type,'fullmatrix')
    C.M = Ac.M*Bc.M;
elseif strcmpi(C.type,'rkmatrix')
    C.U = Ac.U;
    C.V = Bc.V*(Bc.U'*Ac.V);
    C = rSVD_rkmatrix(C,max([Ac.eps Bc.eps]));
end
C.eps = max([Ac.eps Bc.eps]);

% C = hmtx_compress(C);

end


function Hout = formatconversion(Hin,Hout)
    
    if strcmpi(Hin.type,Hout.type)
        Hout = Hin;
    else
        if strcmpi(Hin.type,'supermatrix') && strcmpi(Hout.type,'fullmatrix')
            Hout = super2full(Hin);
        elseif strcmpi(Hin.type,'supermatrix') && strcmpi(Hout.type,'rkmatrix')
            Hout = super2rkmatrix(Hin);
        elseif strcmpi(Hin.type,'rkmatrix') && strcmpi(Hout.type,'fullmatrix')
            Hout = rkmatrix2full(Hin);
        elseif strcmpi(Hin.type,'fullmatrix') && strcmpi(Hout.type,'rkmatrix')
            Hout = full2rkmatrix(Hin);
        elseif strcmpi(Hin.type,'rkmatrix') && strcmpi(Hout.type,'supermatrix')
            Hout = rkmatrix2super(Hin,Hout);
        elseif strcmpi(Hin.type,'fullmatrix') && strcmpi(Hout.type,'supermatrix')
            Hout = full2super(Hin,Hout);
        end
    end
    
end




