% Zachary Chiang, Buenrostro Lab, Harvard University
% code for lamin-ATAC neighborhood analysis in Fig. 3

%% set parameters

home_dir = 'Y:\users\Zack\exigs_code';
table_file = 'TableS1_control_reads.txt';

%% set up environment

cd(home_dir)
addpath(genpath('scripts'));
tic

%% load table

full_table = readtable(sprintf('%s/tables/%s',home_dir,table_file));
pos_col = find(string(full_table.Properties.VariableNames) == 'pos');

disp(sprintf('%s: loaded table',sec2time(toc)))

%% calculate neighborhood correlations

% analysis parameters

dist_thresh = 2;
neighbor_thresh = 5; % minimum number of neighbors

% run for single cell

fov = 0;
cell = 1;

% run on all cells

%cell_list = unique([full_table.fov_idx full_table.cell_idx],'rows');

%for i=1:size(cell_list,1)
    
%fov = cell_list(i,1);
%cell = cell_list(i,2);

disp(sprintf('%s: running fov %d, cell %d',sec2time(toc),fov,cell))

cell_table = full_table(full_table.fov_idx == fov & full_table.cell_idx == cell & full_table.in_nuc==1 & full_table.high_conf_cluster>0,:);

% calculate all distances

dist_mat = pdist2([cell_table.x_um cell_table.y_um cell_table.z_um],[cell_table.x_um cell_table.y_um cell_table.z_um])./exp_factor;

% loop through all ExIGS reads

atac_corr = nan(size(cell_table,1),1);

for i=1:size(cell_table,1);

    % find neighbors and check if enough
    
    neighbor_sel = dist_mat(:,i) < dist_thresh;
    
    if sum(neighbor_sel)>neighbor_thresh
        
        % get neighbor distance to lamin
        
        neighbor_dist = cell_table.dist_to_lamin(neighbor_sel);
        neighbor_dist(neighbor_dist>0.5) = 0.5;
        
        % get neighbor normalized ATAC
        
        neighbor_atac = log2(atac_50kb.Var7(cell_table.bins_50kb(neighbor_sel))./mean(atac_50kb.Var7));

        sel = ~isinf(neighbor_atac);
        atac_corr(i) = corr(neighbor_dist(sel),neighbor_atac(sel),'Rows','pairwise');
        
    end
end

% add to table

cell_table.atac_corr_2um = atac_corr;

%end % to run all cells
