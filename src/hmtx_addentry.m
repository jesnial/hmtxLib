function H = hmtx_addentry(H,ir,jc,x)

% HMTX_ADDENTRY add entries in the specified positions
%   M(ir,jc) = M(ir,jc)+x
%
% USE:
% [H,b] = hmtx_addentry(H,ir,jc,x)
%
% INPUTS:
% 'H': H-matrix
% 'ir': row indices
% 'jc': column indices
% 'x': entry values
%
% OUTPUTS:
% 'H': modified H-matrix
%
% NOTE:
%
% VERSION:
% Date: 25.04.2013
% Copyright(C) 2013: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 26.04.2013: added check on indices before processing
% 26.04.2013: removed already processed rows
% 20.11.2013: fixed bug due to old code leftovers

H = addentry(H,int32(ir),int32(jc),x);

end

function [H,ir,jc,x] = addentry(H,ir,jc,x)

if strcmpi(H.type,'supermatrix')
    % recursive call
    [H.M{1,1},ir,jc,x] = addentry(H.M{1,1},ir,jc,x);
    [H.M{1,2},ir,jc,x] = addentry(H.M{1,2},ir,jc,x);
    [H.M{2,1},ir,jc,x] = addentry(H.M{2,1},ir,jc,x);
    [H.M{2,2},ir,jc,x] = addentry(H.M{2,2},ir,jc,x);
elseif strcmpi(H.type,'fullmatrix')
    if ~isempty(intersect(H.irow,ir(:))) && ~isempty(intersect(H.jcol,jc(:)))
        % expand local indices (reversed to be compliant with linear indexing)
        [jcol,irow] = meshgrid(H.jcol,H.irow);
        % find local matches
        [tf,loc] = ismember([irow(:), jcol(:)],[ir(:), jc(:)],'rows');
        iloc = loc(loc ~= 0);
        % update matrix
        % H.M(tf) = alpha(iloc).*H.M(tf)+x(iloc);
        H.M(tf) = H.M(tf)+x(iloc);
        % remove processed data
        ir(iloc) = [];
        jc(iloc) = [];
        x(iloc) = [];
    end
elseif strcmpi(H.type,'rkmatrix')
    if ~isempty(intersect(H.irow,ir(:))) && ~isempty(intersect(H.jcol,jc(:)))
        % expand local indices (reversed to be compliant with linear indexing)
        [jcol,irow] = meshgrid(H.jcol,H.irow);
        % find local matches
        [tf,loc] = ismember([irow(:), jcol(:)],[ir(:), jc(:)],'rows');
        iloc = loc(loc ~= 0);
        % get subscripts
        [I,J] = ind2sub([H.nrow,H.ncol],find(tf));
        % get dimension
        nel = length(I);
        % add columns to U
        addU = zeros(H.nrow,nel);
        addU(sub2ind([H.nrow,nel],I,(1:nel).')) = x(iloc);
        H.U = [H.U addU];
        % add columns to V
        addV = zeros(H.ncol,nel);
        addV(sub2ind([H.ncol,nel],J,(1:nel).')) = 1;
        H.V = [H.V addV];
        % remove processed data
        ir(iloc) = [];
        jc(iloc) = [];
        x(iloc) = [];
    end
end
end
