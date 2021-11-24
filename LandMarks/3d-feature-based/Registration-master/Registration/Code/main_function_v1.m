% XP debug
% follow the https://github.com/dgoodwin208/Registration/wiki
% First you should download the VLFeat toolbox 

% clear;
addpath(genpath('isoPTVpackage/'));
addpath(genpath('E:\xp\prepocessing'));
basepth = 'E:/data_prj/dir_dataset/DIR_files/';

TRE_e=zeros(10,1);

for idx = 1:1
    pts_struct = DIR_get_all_points_for_the_case(idx, basepth);
    [volmov, spc] = read_DIR_volume_4dCT(idx, 0, basepth);
    volmov = double(volmov);
    pts_mov = pts_struct.extreme.b;
    [volfix, spc] = read_DIR_volume_4dCT(idx, 5, basepth);
    volfix = double(volfix);
    pts_fix = pts_struct.extreme.e;
    init_size=size(volmov);
    volmov=img_thr(volmov,80,900,1);
    volfix=img_thr(volfix,80,900,1);
    min_max1 = [ min(pts_mov, [], 1)', max(pts_mov, [], 1)'];
    min_max2 = [ min(pts_fix, [], 1)', max(pts_fix, [], 1)'];
    min_max = [ min(min_max1(:, 1), min_max2(:, 1)), max(min_max1(:, 2), min_max2(:, 2))];
    d = [10, 10, 5];
    crop_v = [ max(1, min_max(1,1) - d(1)), min(size(volmov, 1), min_max(1,2) + d(1)); ... 
               max(1, min_max(2,1) - d(2)), min(size(volmov, 2), min_max(2,2) + d(2)); ... 
               max(1, min_max(3,1) - d(3)), min(size(volmov, 3), min_max(3,2) + d(3));];
    volmov = crop_data(volmov, crop_v);
    volfix = crop_data(volfix, crop_v);
    volfilename = fullfile('C:\Users\xp\Documents\GitHub\MCMCregistration\preprocessing\Data\',sprintf('case%d_volmov.tif',idx));
    save3DTif(volmov,volfilename);
    fixfilename = fullfile('C:\Users\xp\Documents\GitHub\MCMCregistration\preprocessing\Data\',sprintf('case%d_volfix.tif',idx));
    save3DTif(volfix,fixfilename);
    sample_name=fullfile(sprintf('case%d_',idx));
    calculateDescriptors_4dct(sample_name,1,1,9);
    calculateDescriptors_4dct(sample_name,2,1,9);
    registerWithDescriptors_correct(2,sample_name); 
    output_TPS_dir='C:\Users\xp\Documents\GitHub\MCMCregistration\preprocessing\output';
    output_TPS_filename = fullfile(output_TPS_dir,[sample_name,'TPS.mat']);
    load(output_TPS_filename);
    [r,c,s]=size(volmov);
    T=zeros(r,c,s,3);
    for n=1:length(out1D_total)
       num=in1D_total(n);
       [x_index,y_index,z_index]=ind2sub([r,c,s],num);
       [x_trans,y_trans,z_trans]=ind2sub([r,c,s],out1D_total(n));
       T(x_index,y_index,z_index,1)=x_trans-x_index;
       T(x_index,y_index,z_index,2)=y_trans-y_index;
       T(x_index,y_index,z_index,3)=z_trans-z_index; 
    end

    [volmov_or, Knots_n_rsz] = uncrop_data(volmov, T, crop_v, init_size);

    landregis=Land_trace(Knots_n_rsz,pts_mov,init_size);

    TRE_e(idx)=Get_TREs(landregis,pts_mov,pts_fix,spc);
end



