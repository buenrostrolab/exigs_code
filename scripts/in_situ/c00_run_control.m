% Author: Zachary Chiang, Buenrostro Lab, Harvard University
%
% Processes in situ sequencing and expansion immunofluorescence images 
% for control skin fibroblasts (s83)
%
% This is an example of a wrapper script used to set sample parameters 
% and runs scripts from the main_steps directory
%
% Each step should be able to re-run independently

%% 1. Segment nuclei from base1 of sequencing

% set up environment

home_dir = 'Y:\users\Zack\xgs';
cd(home_dir)
addpath(genpath('scripts'))

% sample parameters

expt = '230530_ctrl_fibroblast';
fov_conv = dlmread(sprintf('ref/fov_conv/%s.txt',expt)); % used to skip fields removed during sequencing
base1_pre = 'Y:\users\IGS\Samples\Sample83_Expt23\Base1\S83_B1'; % everything besides the _F001.ims part of the image name
dapi_thresh = 40;
skip_fovs = [];

% loop through all fields

for fov=0:1:40
    
    if sum(skip_fovs==fov)>0
        continue
    end
    
    adj_fov = fov_conv(fov_conv(:,3)==fov,1); % finds the original field number during base1
    base1_file = sprintf('%s_F%03d.ims',base1_pre,adj_fov);
    c01_seg_base1(home_dir,expt,fov,base1_file,dapi_thresh);
    
end

%% 2. Offset fields from bases 2-n and save cropped nuclei stacks

% set up environment

home_dir = 'Y:\users\Zack\xgs';
cd(home_dir)
addpath(genpath('scripts'))

% sample parameters

expt = '230530_ctrl_fibroblast';
fov_conv = dlmread(sprintf('ref/fov_conv/%s.txt',expt));
file_paths = readtable(sprintf('ref/paths/%s.txt',expt),'Delimiter',',');
base1_pre = 'Y:\users\IGS\Samples\Sample83_Expt23\Base1\S83_B1'; % everything besides the _F001.ims part of the image name
skip_fovs = [];

% loop through all fields and bases (easy to parallelize w/ multiple matlab instances)

for fov=0:1:40
    for base=2:21
    
        if sum(skip_fovs==fov)>0
            continue
        end

        base1_fov = fov_conv(fov_conv(:,3)==fov,1);
        base1_file = sprintf('%s_F%03d.ims',base1_pre,base1_fov);

        % this adjusts for missing fields, disrupted bases, etc.
        
        seq_pre = string(file_paths.img_name(file_paths.base == base));
        if base<=2
            adj_fov = fov_conv(fov_conv(:,3)==fov,1);
            seq_file = sprintf('%s_F%03d.ims',seq_pre,adj_fov);
        elseif base>=3 & base <= 6
            adj_fov = fov_conv(fov_conv(:,3)==fov,2);
            seq_file = sprintf('%s_F%02d.ims',seq_pre,adj_fov);
        else
            seq_file = sprintf('%s_F%02d.ims',seq_pre,fov);
        end
      
        c02_offset_and_crop(home_dir,expt,fov,base,base1_file,seq_file);
    
    end  
end

%% 3. Registration of nuclei stacks

% step 3 (registration) is far faster in on an HPC cluster, see "c00_run_control_reg.m"
% cropped images and nuclei bounds must be transferred to the cluster before registration
% registered images and correlations must be transerred back from the cluster after registration

%% 4a. Peak calling of 3D amplicon locations (originally forgot to number this step, hence 4a and 4b)

% set up environment

home_dir = 'Y:\users\Zack\xgs';
cd(home_dir)
addpath(genpath('scripts'))

% sample parameters

expt = '230530_ctrl_fibroblast';
dapi_thresh = 40;
base1_peak_thresh = [25 25 15 25]; % base 1 (set manually per sample)
base2_peak_thresh = [80 75 80 80]; % base 2 (set manually per sample)
recall_peaks = 1;
voxel_size = [.15 .15 .3];
skip_cells = [];

% loop through fields and cells

for fov=0:1:40
    
    fov_dir = sprintf('processed/%s/fov%03d',expt,fov);
    all_bounds = dlmread(sprintf('%s/cell_bounds.txt',fov_dir));
    
    for cell=1:size(all_bounds,1)
        
        if any(skip_cells(:,1)==fov & skip_cells(:,2)==cell)
            continue
        end
        
        c04_call_peaks(home_dir,expt,fov,cell,dapi_thresh,base1_peak_thresh,recall_peaks,voxel_size)
        c04_call_peaks_base2(home_dir,expt,fov,cell,dapi_thresh,base2_peak_thresh,recall_peaks,voxel_size)
    end  
end

%% 4b. Correction of SBS signal

% set up environment

home_dir = 'Y:\users\Zack\xgs';
cd(home_dir)
addpath(genpath('scripts'))

% sample parameters

expt = '230530_ctrl_fibroblast';
bases = [1:6 7:9 13:18 19:21]';
seq_peak_thresh = [25 25 15 25]; % copy from 4a
voxel_size = [.15 .15 3];
skip_cells = [];

% loop through all fields and cells

for fov=0:1:40
    
    fov_dir = sprintf('processed/%s/fov%03d',expt,fov);
    all_bounds = dlmread(sprintf('%s/cell_bounds.txt',fov_dir));
    
    for cell=1:size(all_bounds,1)
    
        cell_dir = mkdirp(sprintf('%s/cell%03d',fov_dir,cell));

        if any(skip_cells(:,1)==fov & skip_cells(:,2)==cell)
            continue
        end

        for base=1:21

            % make sure registered images exist
            
            if ~exist(sprintf('%s/reg/base%02d_ch4.tif',cell_dir,base))
                continue
            end

            if base>1
                prev_base = bases(find(bases==base)-1);
                c04_correct_sbs(home_dir,expt,fov,cell,base,prev_base,seq_peak_thresh,voxel_size)
            else
                c04_correct_mixing(home_dir,expt,fov,cell,base,0,seq_peak_thresh,voxel_size)
            end
        end
    end
end

%% 5. Matching of in situ and ex situ UMIs (see "ex_situ" directories)

% set up environment

home_dir = 'Y:\users\Zack\xgs';
cd(home_dir)
addpath(genpath('scripts'))

% sample parameters

home_dir = 'Y:\users\Zack\xgs';
expt = '230530_ctrl_fibroblast';
sel_bases = [1:9 13:21];
seq_peak_thresh = [25 25 15 25; 80 75 80 80]; % [base1; base2], can take from step 4a or adjust here
voxel_size = [.15 .15 .3];
umi_pre = sprintf('%s/ex_situ/%s/230530_sample83_#_next2_nova1.umis_all.repeats.txt',home_dir,expt); % "#" in place of ex situ sample number
num_umi_samples = 4;

skip_table = readtable(sprintf('ref/skip_idx/%s.txt',expt)); % used to skip particular bases for certain cells
skip_cells = [];

% loop through all fields and cells

for fov=0:1:40
    
    fov_dir = sprintf('processed/%s/fov%03d',expt,fov);
    cell_bounds = dlmread(sprintf('%s/cell_bounds.txt',fov_dir));
    num_cells = size(cell_bounds,1);
    
    for cell=1:num_cells;
        
        if any(skip_cells(:,1)==fov & skip_cells(:,2)==cell)
            continue
        end
        
        cell_dir = mkdirp(sprintf('%s/cell%03d',fov_dir,cell));
        skips = str2num(skip_table.skips{skip_table.fov == fov & skip_table.cell == cell});
        sel_bases = sel_bases(~ismember(sel_bases,skips));
        
        c05_match_umis(home_dir,expt,fov,cell,seq_peak_thresh,sel_bases,voxel_size,umi_pre,num_umi_samples)

    end
end

%% 6b. Quantify expansion factor
% before running this step, you must first find the affine transform between pre- and post-expansion montages
% this is done once per sample and is easiest interactively which is why it is not called from here
% see c06a_register_pre_post_montages.m for an example

% set up environment

home_dir = 'Y:\users\Zack\xgs';
cd(home_dir)
addpath(genpath('scripts'))

% sample parameters

expt = '230530_sample83_v3';
post_offset = [0 0];
base1_pre = 'Y:\users\IGS\Samples\Sample83_Expt23\Base1\S83_B1';
post_exp_img = '';
pre_exp_img = 'Y:\users\IGS\Experiments\Expt23_NP_HP1HP2_linking\10X_PreEX_montage\10XPreEx_montage_Expt23_FusionStitcher.ims';
coords_file = 'Y:\users\IGS\Samples\Sample83_Expt23\Base1\230530_sample83_v3.coords.txt';
skip_fovs = [];

% loop through all fields and cells

for fov=0:1:40
    
    fov_dir = sprintf('processed/%s/fov%03d',expt,fov);
    cell_bounds = dlmread(sprintf('%s/cell_bounds.txt',fov_dir));
    num_cells = size(cell_bounds,1);
    
    for cell=1:num_cells;
        
        cell_dir = sprintf('%s/cell%03d',fov_dir,cell);
        if any(skip_fovs(:,1)==fov & skip_fovs(:,2)==cell)
            continue
        end

        c09_cell_exp_factor_v3(home_dir,expt,fov,cell,post_offset,base1_pre,post_exp_img,pre_exp_img,coords_file)
        
    end
end

%% 7. Offset expasion IF and save cropped nuclei stacks (like 2 but for IF)

% set up environment

home_dir = 'Y:\users\Zack\xgs';
cd(home_dir)
addpath(genpath('scripts'))

% sample parameters

expt = '230530_ctrl_fibroblast';
fov_conv = dlmread(sprintf('ref/fov_conv/%s.txt',expt));
base1_pre = 'Y:\users\IGS\Samples\Sample83_Expt23\Base1\S83_B1';
if_pre = 'Y:\users\IGS\Samples\Sample83_Expt23\LMNA_IF\S83_LMNA';
skip_fovs = [];

% loop through all fields

for fov=0:1:40
    
    if any(skip_fovs==fov)
        continue
    end
    
    adj_fov = fov_conv(fov_conv(:,3)==fov,1);
    base1_file = sprintf('%s_F%03d.ims',base1_pre,adj_fov);
    if_file = sprintf('%s_F%03d.ims',if_pre,adj_fov);

    c07_offset_if(home_dir,expt,fov,base1_file,if_file)
    
end

%% 8. Segmentation of peripheral and internal lamin

% set up environment

home_dir = 'Y:\users\Zack\xgs';
cd(home_dir)
addpath(genpath('scripts'))

% sample parameters

expt = '230530_ctrl_fibroblast';
exp_factor = 4.7; % calculate from step 6
resegment = 1; % set to 0 if thresholds haven't changed
skip_cells = [];
if_thresh = readtable(sprintf('ref/if_thresh_list/%s.txt',expt)); % occasionally used for manual adjustment

% loop through all fields and cells

for fov=0:1:40
    
    fov_dir = sprintf('processed/%s/fov%03d',expt,fov);
    cell_bounds = dlmread(sprintf('%s/cell_bounds.txt',fov_dir));
    num_cells = size(cell_bounds,1);
    
    for cell=1:num_cells;
        
        if any(skip_cells(:,1)==fov & skip_cells(:,2)==cell)
            continue
        end
        
        seg_thresh = if_thresh{if_thresh.fov == fov & if_thresh.cell == cell,3:4};
        c08_seg_lamin(home_dir,expt,fov,cell,seg_thresh(1),seg_thresh(2),exp_factor,resegment)

    end
end

%% 9. Cluster reads into homologous chromsomes

% sample environment

home_dir = 'Y:\users\Zack\xgs';
cd(home_dir)
addpath(genpath('scripts'))

% sample parameters

expt = '230530_ctrl_fibroblast';
exp_factor = 4.7; % calculate from step 6
table_file = 'match_tables/lamin_table_v9.txt';
skip_cells = [];

% loop through all fields and cells

for fov=0:40
    
    fov_dir = sprintf('processed/%s/fov%03d',expt,fov);
    cell_bounds = dlmread(sprintf('%s/cell_bounds.txt',fov_dir));
    num_cells = size(cell_bounds,1);
    
    for cell=1:num_cells;
        
        if any(skip_cells(:,1)==fov & skip_cells(:,2)==cell)
            continue
        end
        
        c08_cluster_chrs_v2(home_dir,expt,fov,cell,exp_factor,table_file)

    end
end



