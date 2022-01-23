function [HT] = cluster_for_hm(N_leaf,leaf,bar,eta)
% create a cluster tree for sorted points according to distance (not
% balanced according to the number of points)
%
% VARIABLES :
% 'Nleaf ' : number of leaves of the cluster tree
% 'leaf' : leaves of the cluster tree
% 'bar' : cluster of points
% 'eta' : admissibility criterion
%


% change variable names

HT=hmatrix(zeros(2));
HT.F=[]; %create matrix
HT.admissible=0; % only for leaf nodes
HT.sz=[sum(leaf(2,:)) sum(leaf(2,:))]; % matrix size

% center and diameters of rows and columns
[Crows_center,Crows_diam] = fun_my_center_and_radius_cluster(bar);
[Ccolumns_center,Ccolumns_diam] = fun_my_center_and_radius_cluster(bar);

% separate clutser into 4
HT=fun_my_create_4sons(N_leaf,N_leaf,leaf,leaf,HT,Crows_center,Crows_diam,Ccolumns_center,Ccolumns_diam,bar,bar,eta);
end 

function [HT] = fun_my_create_4sons(Nleaf_rows,Nleaf_columns,leaf_rows,leaf_columns,HT,Crows_center,Crows_diam,Ccolumns_center,Ccolumns_diam,bar_rows,bar_columns,eta)

     % separate cluster into 4 blocks : create 4 sons in the cluster tree
     %
     % VARIABLES : 
     % 'Nleaf_rows': number of leaves in the rows
     % 'Nleaf_columns': number of leaves in the columns
% 
%      Nleaf_rows_left=floor(Nleaf_rows/2); % change separation method here
%      Nleaf_rows_right=ceil(Nleaf_rows/2);
%      Nleaf_columns_left=floor(Nleaf_columns/2);
%      Nleaf_columns_right=ceil(Nleaf_columns/2);

    % get splitting location
    middle_rows = abs(min(bar(1,:))+Crows_diam/2);
    middle_columns = abs(min(bar(1,:))+Ccolumns_diam/2);
    
    % get splitting points
    Nleaf_rows_left = sum(bar<middle_rows);
    Nleaf_rows_right = Nleaf_rows - Nleaf_rows_left;
    Nleaf_columns_left = sum(bar<middle_columns);
    Nleaf_columns_right = Nleaf_columns - Nleaf_columns_left;
    
    
     % get block sizes
     ssize_rows_left=sum(leaf_rows(2,1:Nleaf_rows_left));
     ssize_rows_right=sum(leaf_rows(2,Nleaf_rows_left+1:end));
     ssize_columns_left=sum(leaf_columns(2,1:Nleaf_columns_left));
     ssize_columns_right=sum(leaf_columns(2,Nleaf_columns_left+1:end));
%

    % get leaves for each block
    leaf_rows_left=leaf_rows(:,1:Nleaf_rows_left);
    leaf_rows_right=leaf_rows(:,Nleaf_rows_left+1:end);
    leaf_columns_left=leaf_columns(:,1:Nleaf_columns_left);
    leaf_columns_right=leaf_columns(:,Nleaf_columns_left+1:end);
%
    % get sub-clusters 
    bar_rows_left=bar_rows(:,1:ssize_rows_left); 
    bar_rows_right=bar_rows(:,ssize_rows_left+1:end);
    bar_columns_left=bar_columns(:,1:ssize_columns_left);
    bar_columns_right=bar_columns(:,ssize_columns_left+1:end);
%
    % get center and radius for each block 
    [Crows_center_left,Crows_diam_left] = fun_my_center_and_radius_cluster(bar_rows_left);
    [Crows_center_right,Crows_diam_right] = fun_my_center_and_radius_cluster(bar_rows_right);
    [Ccolumns_center_left,Ccolumns_diam_left] = fun_my_center_and_radius_cluster(bar_columns_left);
    [Ccolumns_center_right,Ccolumns_diam_right] = fun_my_center_and_radius_cluster(bar_columns_right);
 %
    % create empty matrix blocks 
    
     HT.A11=hmatrix(zeros(2));HT.A11.F=[];
     HT.A11.sz=[ssize_rows_left,ssize_columns_left];
     HT.A11.admissible=0;
 %
     HT.A12=hmatrix(zeros(2));HT.A12.F=[];
     HT.A12.sz=[ssize_rows_left,ssize_columns_right];  
     HT.A12.admissible=0;
 %
     HT.A21=hmatrix(zeros(2));HT.A21.F=[];
     HT.A21.sz=[ssize_rows_right,ssize_columns_left];
     HT.A21.admissible=0; 
 % 
     HT.A22=hmatrix(zeros(2));HT.A22.F=[];
     HT.A22.sz=[ssize_rows_right,ssize_columns_right];
     HT.A22.admissible=0;
 %     
 
    % recursive part
    % call function again for admissible blocks
     if fun_check_adm(Crows_center_left,Crows_diam_left,Ccolumns_center_left,...
                    Ccolumns_diam_left,eta) % check for admissibility
        HT.A11.admissible=1; %admissible block
     elseif  Nleaf_rows_left>1 && Nleaf_columns_left>1 % non admissible block : if block size > Nmin (here Nmin =1, change this)
            [HT.A11] = fun_my_create_4sons(Nleaf_rows_left,Nleaf_columns_left,...
           leaf_rows_left,leaf_columns_left,HT.A11,Crows_center,Crows_diam,...
            Ccolumns_center,Ccolumns_diam,bar_rows_left,bar_columns_left,eta); % create sub blocks
        % else : it's a full block
     end
     

     if fun_check_adm(Crows_center_left,Crows_diam_left,Ccolumns_center_right,Ccolumns_diam_right,eta) 
        HT.A12.admissible=1;
     elseif Nleaf_rows_left>1 && Nleaf_columns_right>1
        [HT.A12] = fun_my_create_4sons(Nleaf_rows_left,Nleaf_columns_right,...
              leaf_rows_left,leaf_columns_right,HT.A12,Crows_center,Crows_diam,...
          Ccolumns_center,Ccolumns_diam,bar_rows_left,bar_columns_right,eta);
     end

     if fun_check_adm(Crows_center_right,Crows_diam_right,Ccolumns_center_left,...
                Ccolumns_diam_left,eta) 
         HT.A21.admissible=1;
     elseif Nleaf_rows_right>1 && Nleaf_columns_left>1
        [HT.A21] = fun_my_create_4sons(Nleaf_rows_right,Nleaf_columns_left,...
                 leaf_rows_right,leaf_columns_left,HT.A21,Crows_center,Crows_diam,...
                Ccolumns_center,Ccolumns_diam,bar_rows_right,bar_columns_left,eta);
     end

     if fun_check_adm(Crows_center_right,Crows_diam_right,Ccolumns_center_right,...
                 Ccolumns_diam_right,eta) 
        HT.A22.admissible=1;
     elseif Nleaf_rows_right>1 && Nleaf_columns_right>1
        [HT.A22] = fun_my_create_4sons(Nleaf_rows_right,Nleaf_columns_right,...
                    leaf_rows_right,leaf_columns_right,HT.A22,Crows_center,Crows_diam,...
                    Ccolumns_center,Ccolumns_diam,bar_rows_right,bar_columns_right,eta);
     end


end

function [C_center,C_diam] = fun_my_center_and_radius_cluster(bar)
%get center (mean) and diameter of the cluster
%
% VARIABLES : 
% 'bar' : cluster of points

 Ndofs=size(bar,2);
    C_center=sum(bar,2)/Ndofs; %center
    C_radii =sqrt((bar(1,:)-C_center(1)).^2+...
                  (bar(2,:)-C_center(2)).^2+...
                  (bar(3,:)-C_center(3)).^2); 
    C_diam=2*max(C_radii); %diameter
end
