% Zachary Chiang, Buenrostro Lab, Harvard University
% code to compare % lamin and total lamin-ATAC correlations

%% set parameters

home_dir = 'Y:\users\Zack\exigs_code';
control_table_file = 'TableS1_control_reads.txt';
progeria_table_file = 'TableS3_progeria_reads.txt';
stats_table_file = 'TableS4_cell_stats.txt';

%% set up environment

cd(home_dir)
addpath(genpath('scripts'));
tic

%% load table

control_table = readtable(sprintf('%s/tables/%s',home_dir,control_table_file));
progeria_table = readtable(sprintf('%s/tables/%s',home_dir,progeria_table_file));
cell_stats_table = readtable(sprintf('%s/tables/%s',home_dir,stats_table_file));

disp(sprintf('%s: loaded tables',sec2time(toc)))

%% calculate total lamin-ATAC corr for control cells

control_stats = cell_stats_table(string(cell_stats_table.expt_name) == 'progeria_ctrl',:);
control_cells = unique(control_stats{:,2:3},'rows');
num_control_cells = size(control_cells,1);
control_corr = nan(num_control_cells,1);

for cell_idx=1:num_control_cells
    
    
    sel_fov = control_cells(cell_idx,1);
    sel_cell = control_cells(cell_idx,2);
    
    cell_table = control_table(control_table.fov_idx == sel_fov & control_table.cell_idx == sel_cell ...
        & control_table.in_nuc & control_table.high_conf_cluster > 0,:);
    
    if size(cell_table,1) == 0
        continue
    end
    
    dist_to_lamin = cell_table.dist_to_lamin;
    dist_to_lamin(dist_to_lamin>0.5) = 0.5;
    atac = cell_table.atac_50kb;
    sel = ~isinf(atac);
    
    control_corr(cell_idx) = corr(dist_to_lamin(sel),atac(sel),'Rows','pairwise');
    
end

%% calculate total lamin-ATAC corr for progeria cells

progeria_stats = cell_stats_table(string(cell_stats_table.expt_name) == 'progeria_p19' ...
    | string(cell_stats_table.expt_name) == 'progeria_p22' | string(cell_stats_table.expt_name) == 'progeria_p25',:);
progeria_cells = unique(progeria_stats(:,1:3),'rows');
num_progeria_cells = size(progeria_cells,1);
progeria_corr = nan(num_progeria_cells,1);

for cell_idx=1:num_progeria_cells
    
    sel_expt = progeria_cells{cell_idx,1};
    sel_fov = progeria_cells{cell_idx,2};
    sel_cell = progeria_cells{cell_idx,3};
    
    cell_table = progeria_table(progeria_table.fov_idx == sel_fov & progeria_table.cell_idx == sel_cell ...
        & progeria_table.in_nuc & progeria_table.high_conf_cluster > 0,:);
    
    if size(cell_table,1) == 0
        continue
    end
    
    dist_to_lamin = cell_table.dist_to_lamin;
    dist_to_lamin(dist_to_lamin>0.5) = 0.5;
    atac = cell_table.atac_50kb;
    sel = ~isinf(atac);
    
    progeria_corr(cell_idx) = corr(dist_to_lamin(sel),atac(sel),'Rows','pairwise');
    
end

%% compare % lamin and lamin-ATAC corr

progeria_samples = ["progeria_p19","progeria_p22","progeria_p25"];

figure;
scatter(control_corr,control_stats.prc_lam); hold on;
for i=1:size(progeria_samples,2)
    sel = progeria_stats.expt_name == progeria_samples(i);
    scatter(progeria_corr(sel),progeria_stats.prc_lam(sel)); hold on;
end

xlabel('total lamin-ATAC corr'); ylabel('% lamin')
xlim([-0.1 0.4]); ylim([0 0.5])

%% compare % internal and lamin-ATAC corr

progeria_samples = ["progeria_p19","progeria_p22","progeria_p25"];

figure;
scatter(control_corr,control_stats.prc_int); hold on;
for i=1:size(progeria_samples,2)
    sel = progeria_stats.expt_name == progeria_samples(i);
    scatter(progeria_corr(sel),progeria_stats.prc_int(sel)); hold on;
end

xlabel('total lamin-ATAC corr'); ylabel('% internal lamin')
xlim([-0.1 0.4]); ylim([0 0.15])

