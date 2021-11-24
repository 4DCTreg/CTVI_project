% follow the https://github.com/dgoodwin208/Registration/wiki
% First you should download the VLFeat toolbox

clear;

addpath(genpath('../../VAMPIRE_FullDatabase_MHA/'));
addpath(genpath('./scr/'));
addpath(genpath('./pTVreg-master/mutils/My/'));
addpath(genpath('./pTVreg-master/ptv'));
addpath(genpath('./Reg_method/'));

% RefVI_type = 'Galligas';
% RefVI_type = 'Xenon';
RefVI_type = 'DTPA-SPECT';
% VI_metric = 'DIR_Jac';
VI_metric = 'DIR_HU';

if strcmp(RefVI_type,'Galligas')
    file_path = 'H:/VAMPIRE_FullDatabase_MHA/Study01_Galligas-PET/';
    landmarks_path = 'H:/CTVI/LandMarks/SD_01_v1/';
    subject_num = 25;
    num_phase = 5;
else
    if strcmp(RefVI_type,'Xenon')
        file_path = 'H:/VAMPIRE_FullDatabase_MHA/Study02_Xenon-CT/';
        landmarks_path = 'H:/CTVI/LandMarks/SD_02/';
        subject_num = 4;
        num_phase = 8;
    else
        if strcmp(RefVI_type,'DTPA-SPECT')
            file_path = 'H:/VAMPIRE_FullDatabase_MHA/Study03_DTPA-SPECT/';
            landmarks_path = 'H:/CTVI/LandMarks/SD_03/';
            subject_num = 21;
            num_phase = 10;
        end
    end
end
TREs = zeros(1, subject_num);
TREs_or = zeros(1, subject_num);
SpR = zeros(1, subject_num);


save_path = 'H:/CTVI/Results/Continuous_figure/SD03/';
ex_to_avg_path = 'H:/CTVI/Results/SD_03_ex_to_avg_sup_multi/';
for num = 6:subject_num
    if num < 10
        CT_path = [file_path,'Subject_0',num2str(num),'/CT/'];
        Ref_path = [file_path,'Subject_0',num2str(num),'/GT/'];
    else
        CT_path = [file_path,'Subject_',num2str(num),'/CT/'];
        Ref_path = [file_path,'Subject_',num2str(num),'/GT/'];
    end
    Ref_Filename = [Ref_path,'VentImage.mha'];
    info = mha_read_header(Ref_Filename);
    GT_img = mha_read_volume(info);
    GT_mask_filename = [Ref_path,'VentMask.mha'];
    info_mask_GT = mha_read_header(GT_mask_filename);
    GT_mask = mha_read_volume(info_mask_GT);
    
    vol_ct = [];
    vol_mask = [];
    if strcmp(RefVI_type,'Xenon')
        for phase = 1:num_phase
            if phase < 10
                CT_Mask_Filename = [CT_path,'PhaseMask_0',num2str(phase),'.mha'];
            else
                CT_Mask_Filename = [CT_path,'PhaseMask_',num2str(phase),'.mha'];
            end
            info_mask = mha_read_header(CT_Mask_Filename);
            vol_mask_cur = mha_read_volume(info_mask);
            vol_mask = cat(4,vol_mask,vol_mask_cur);
        end
    end
    if strcmp(RefVI_type,'Xenon')
        avg_img_orig = [];
    else
        avg_img_name = [CT_path,'AverageImage.mha'];
        info = mha_read_header(avg_img_name);
        avg_img_orig = mha_read_volume(info);
        avg_mask_name = [CT_path,'AverageMask.mha'];
        info_mask = mha_read_header(avg_mask_name);
        avg_mask_orig = mha_read_volume(info_mask);
    end
    coronal = [];
    for phase = 1:num_phase
        if phase < 10
            CT_Filename = [CT_path,'PhaseImage_0',num2str(phase),'.mha'];
        else
            CT_Filename = [CT_path,'PhaseImage_',num2str(phase),'.mha'];
        end
        if num == 20 && phase==1
            info = mha_read_header(CT_Filename);
            vol_ct_cur = mha_read_volume(info);
            vol_20_01 = vol_ct_cur;
            vol_ct_cur = imresize3(vol_ct_cur, [512,512,84], 'nearest');
        else
            info = mha_read_header(CT_Filename);
            vol_ct_cur = mha_read_volume(info);
        end
        vol_ct = cat(4,vol_ct,vol_ct_cur);
    end
    
    vol_in_orig = double(vol_ct(:,:,:,1));
    vol_ex_orig = double(vol_ct(:,:,:,6));
    init_size = size(vol_in_orig);
    pts_struct = DIR_get_landmarks_for_the_case_03(num, landmarks_path);
    pts_in = pts_struct.extreme.in;
    pts_ex = pts_struct.extreme.ex;
    
    units = info_mask.PixelDimensions;
    
    % homrf reg
    parameter.resize = 1;
    parameter.useTop = 0;
    parameter.metric = 'MIND';
    spc_orig = units;
    [pt_errs_phys_or, TRE_phys_or, TREstd_phys_or] = DIR_or(pts_in, pts_ex, spc_orig, []);
    TREs_or(num) = TRE_phys_or;
    
    if parameter.resize
        bszv = size(vol_in_orig);
        spc_tmp = [1, 1, 1];
        vol_ex = volresize(vol_ex_orig, round(bszv .* units .* spc_tmp), 1);
        vol_in = volresize(vol_in_orig, round(bszv .* units .* spc_tmp), 1);
        avg_mask = imresize3(avg_mask_orig, round(bszv .* units .* spc_tmp), 'nearest');
        spc = [1,1,1] ./ spc_tmp;
    end
    d = [10, 10, 5];
    vol_rsz = size(vol_ex);
    crop_v = crop_mask(avg_mask,size(avg_mask),d);
    sup_real_sz = [crop_v(1,2)-crop_v(1,1)+1,crop_v(2,2)-crop_v(2,1)+1,crop_v(3,2)-crop_v(3,1)+1];
    vol_ex = crop_data(vol_ex,crop_v);
    vol_in = crop_data(vol_in,crop_v);
    vol_ex_ROI = vol_ex;
    vol_in_ROI = vol_in;
    n_seeds = 2000;
    distance = round(nthroot(prod(sup_real_sz)/n_seeds,3));
    h_distance = floor(distance/2);
    [new_crop,min_h_distance] = mulit_sup_corp(crop_v,vol_rsz,h_distance);
    % offset_range
    [offset_Y,offset_X,offset_Z]= meshgrid(-min_h_distance(2):min_h_distance(2),-min_h_distance(1):min_h_distance(1),...
        -min_h_distance(3):min_h_distance(3));
    
    
    %     show_overlay(avg_img_orig, Glandmarks_SIFTT_img, GT_mask, 30, 'autumn')
    vol_in_nonor = vol_in;
    vol_ex_nonor = vol_ex;
    vol_in = img_thr(vol_in, -1024, 0, 1);
    vol_ex = img_thr(vol_ex, -1024, 0, 1);
    parameter.spc = spc;
    Knots_pre_level = zeros(size(vol_ex,1),size(vol_ex,2),size(vol_ex,3),3);
    for CT_phase = 1:num_phase
        if CT_phase == 3
            continue;
        else
            im_index = num_phase - CT_phase + 1;
            vol_mov = vol_ct(:,:,:,im_index);

            if parameter.resize == 1
                vol_mov = volresize(vol_mov, round(bszv .* units .* spc_tmp), 1);
            end
            vol_mov = (crop_data(vol_mov, crop_v));
            vol_mov = img_thr(vol_mov, -1024, 0, 1);
            [Knots_n_tot, Knots_n] = reg_homrf_multi_level_continuous_SD_03(parameter,vol_mov,vol_ex,...
                init_size,pts_ex,pts_in,bszv,crop_v,Knots_pre_level);
            Knots_pre_level = Knots_n_tot;
            savename = [save_path,'homrf_subject_',num2str(num),num2str(im_index),'.mat'];
            save(savename,'Knots_n_tot');
        end
    end
    

end
