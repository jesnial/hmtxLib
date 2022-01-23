function C = hmtx_muladd(C,A,B)

% HMTX_MULADD multiplies two H-matrices and sum the result to another
% H-matrix
%
% USE:
% C = hmtx_muladd(C,A,B)
%
% INPUTS:
% 'C': H-matrix to add
% 'A': H-matrix
% 'B': H-matrix
%
% OUTPUTS:
% 'C': H-matrix such that C = C+A*B
%
% NOTE:
% The result of the product has the most convenient matrix format
%
% VERSION:
% Date: 04.02.2014
% Copyright(C) 2014: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 06.02.2014: added compatibility check
% 06.02.2014: refurbished routine

% check on block cluster tree
if ~isequal(A.jcol,B.irow)
    error('hmtx_muladd: H-matrices not compatible\n')
end

if ~isempty(C) && strcmpi(C.type,'supermatrix') && strcmpi(A.type,'supermatrix') && strcmpi(B.type,'supermatrix')
    for i = 1:2
        for j = 1:2
            for k = 1:2
                C.M{i,j} = hmtx_muladd(C.M{i,j},A.M{i,k},B.M{k,j});
            end
        end
    end
else
    if strcmpi(A.type,'supermatrix') && strcmpi(B.type,'supermatrix')
        AxB = hmtx_create('supermatrix',A.irow,B.jcol);
        for i = 1:2
            for j = 1:2
                for k = 1:2
                    AxB.M{i,j} = hmtx_muladd(AxB.M{i,j},A.M{i,k},B.M{k,j});
                end
            end
        end
    elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'fullmatrix')
        % create A*B H-matrix
        AxB = hmtx_create('fullmatrix',A.irow,B.jcol);
        AxB.M = A.M*B.M;
        AxB.eps = max([A.eps,B.eps]);
    elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'rkmatrix')
        % create A*B H-matrix
        AxB = hmtx_create('rkmatrix',A.irow,B.jcol);
        if A.k ~= 0 && B.k ~= 0
            AxB.U = A.U;
            AxB.V = B.V*(B.U'*A.V);
            AxB.eps = max([A.eps,B.eps]);
            AxB.k = A.k;
        end
        % recompression
        %         AxB = rSVD_rkmatrix(AxB,AxB.eps);
    elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'fullmatrix')
        % create A*B H-matrix
        AxB = hmtx_create('rkmatrix',A.irow,B.jcol);
        if A.k ~= 0
            AxB.U = A.U;
            AxB.V = B.M'*A.V;
            AxB.eps = max([A.eps,B.eps]);
            AxB.k = A.k;
        end
        % recompression
        %         AxB = rSVD_rkmatrix(AxB,AxB.eps);
    elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'rkmatrix')
        % create A*B H-matrix
        AxB = hmtx_create('rkmatrix',A.irow,B.jcol);
        if B.k ~= 0
            AxB.U = A.M*B.U;
            AxB.V = B.V;
            AxB.eps = max([A.eps,B.eps]);
            AxB.k = B.k;
        end
        % recompression
        %         AxB = rSVD_rkmatrix(AxB,AxB.eps);
    elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'rkmatrix')
        if B.k ~= 0 && (isempty(C) || ~strcmpi(C.type,'supermatrix'))
            % convert A to fullmatrix
            A = super2full(A);
            AxB = hmtx_muladd([],A,B);
        elseif B.k ~= 0
            B = rkmatrix2super(B,A.M{1,1}.jcol,C.M{1,1}.jcol);
            AxB = hmtx_create('supermatrix',A.irow,B.jcol);
            for i = 1:2
                for j = 1:2
                    for k = 1:2
                        AxB.M{i,j} = hmtx_muladd(AxB.M{i,j},A.M{i,k},B.M{k,j});
                    end
                end
            end
        else
            AxB = hmtx_create('rkmatrix',A.irow,B.jcol);
        end
    elseif strcmpi(A.type,'supermatrix') && strcmpi(B.type,'fullmatrix')
        if isempty(C) || ~strcmpi(C.type,'supermatrix')
            % convert A to rkmatrix
            A = super2full(A);
            AxB = hmtx_muladd([],A,B);
        else
            B = full2super(B,A.M{1,1}.jcol,C.M{1,1}.jcol);
            AxB = hmtx_create('supermatrix',A.irow,B.jcol);
            for i = 1:2
                for j = 1:2
                    for k = 1:2
                        AxB.M{i,j} = hmtx_muladd(AxB.M{i,j},A.M{i,k},B.M{k,j});
                    end
                end
            end
        end
    elseif strcmpi(A.type,'rkmatrix') && strcmpi(B.type,'supermatrix')
        if A.k ~= 0 && (isempty(C) || ~strcmpi(C.type,'supermatrix'))
            % convert B to fullmatrix
            B = super2full(B);
            AxB = hmtx_muladd([],A,B);
            %             % convert B to rkmatrix
            %             B = super2rkmatrix(B);
            %             AxB = hmtx_muladd2([],A,B);
        elseif A.k ~= 0
            A = rkmatrix2super(A,C.M{1,1}.irow,B.M{1,1}.irow);
            AxB = hmtx_create('supermatrix',A.irow,B.jcol);
            for i = 1:2
                for j = 1:2
                    for k = 1:2
                        AxB.M{i,j} = hmtx_muladd(AxB.M{i,j},A.M{i,k},B.M{k,j});
                    end
                end
            end
        else
            AxB = hmtx_create('rkmatrix',A.irow,B.jcol);
        end
    elseif strcmpi(A.type,'fullmatrix') && strcmpi(B.type,'supermatrix')
        if isempty(C) || ~strcmpi(C.type,'supermatrix')
            % convert A to rkmatrix
            B = super2full(B);
            AxB = hmtx_muladd([],A,B);
        else
            A = full2super(A,C.M{1,1}.irow,B.M{1,1}.irow);
            AxB = hmtx_create('supermatrix',A.irow,B.jcol);
            for i = 1:2
                for j = 1:2
                    for k = 1:2
                        AxB.M{i,j} = hmtx_muladd(AxB.M{i,j},A.M{i,k},B.M{k,j});
                    end
                end
            end
        end
    end
    % sum with C
    if isempty(C)
        C = AxB;
    else
        C = hmtx_add(C,AxB,1,1);
    end
end

