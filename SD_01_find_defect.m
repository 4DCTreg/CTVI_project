% follow the https://github.com/dgoodwin208/Registration/wiki
% First you should download the VLFeat toolbox

clear;

addpath(genpath('../../VAMPIRE_FullDatabase_MHA/'));
addpath(genpath('./scr/'));
addpath(genpath('G:/SIFT3D/SIFT3D 1.4.5/'));
addpath(genpath('./Reg_method/'))
addpath(genpath('./pTVreg-master/mutils/My/'));
addpath(genpath('./pTVreg-master/ptv'));
ex_to_avg_path = 'H:/CTVI/Results/SD_01_ex_to_avg/';
in_to_ex_path_homrf = 'H:/CTVI/Results/SD_01_in_to_ex_homrf/';
RefVI_type = 'Galligas';
% RefVI_type = 'Xenona';
% RefVI_type = 'DTPA-SPECT';
% VI_metric = 'DIR_HU';
VI_metric = 'DIR_Jac';

if strcmp(RefVI_type,'Galligas')
    file_path = 'H:/VAMPIRE_FullDatabase_MHA/Study01_Galligas-PET/';
    landmarks_path = 'H:/CTVI/LandMarks/SD_01_v1/';
    subject_num = 25;
    num_phase = 5;
else
    if strcmp(RefVI_type,'Xenon')
        file_path = 'H:/VAMPIRE_FullDatabase_MHA/Study02_Xenon-CT/';
        subject_num = 4;
        num_phase = 8;
    else
        if strcmp(RefVI_type,'DTPA-SPECT')
            file_path = 'H:/VAMPIRE_FullDatabase_MHA/Study03_DTPA-SPECT/';
            subject_num = 21;
            num_phase = 10;
        end
    end
end
TREs = zeros(1,25);
TREs_or = zeros(1,25);
SpR = zeros(1, 25);
for num = 25:subject_num
    if num < 10
        CT_path = [file_path,'Subject_0',num2str(num),'/CT/'];
        Ref_path = [file_path,'Subject_0',num2str(num),'/GT/'];
    else
        CT_path = [file_path,'Subject_',num2str(num),'/CT/'];
        Ref_path = [file_path,'Subject_',num2str(num),'/GT/'];
    end
    Ref_Filename = [Ref_path,'VentImage.mha'];
    info_GT = mha_read_header(Ref_Filename);
    GT_img = mha_read_volume(info_GT);
    GT_mask_filename = [Ref_path,'VentMask.mha'];
    info = mha_read_header(GT_mask_filename);
    GT_mask = mha_read_volume(info);
    
    vol_ct = [];
    if strcmp(RefVI_type,'Xenon')
        avg_img_orig = [];
    else
        avg_img_name = [CT_path,'AverageImage.mha'];
        info = mha_read_header(avg_img_name);
        avg_img_orig = mha_read_volume(info);
        avg_mask_name = [CT_path,'AverageMask.mha'];
        info_mask = mha_read_header(avg_mask_name);
        avg_mask = mha_read_volume(info_mask);
    end
    for phase = 1:num_phase
        if phase < 10
            CT_Filename = [CT_path,'PhaseImage_0',num2str(phase),'.mha'];
        else
            CT_Filename = [CT_path,'PhaseImage_',num2str(phase),'.mha'];
        end
        info = mha_read_header(CT_Filename);
        vol_ct_cur = mha_read_volume(info);
        vol_ct = cat(4,vol_ct,vol_ct_cur);
    end
    sz_GT = info_GT.Dimensions;
    sz_ct = info.Dimensions;
    index_slice = 275;
    figure 
    show_overlay_sagittal(avg_img_orig, GT_img, info_GT, GT_mask, index_slice,'hot')
    title_name = ['overlay',num2str(num)];
    title(title_name);
%     d = [10, 10, 5];
%     crop_v = crop_mask(avg_mask,info_mask,d);
%     vol_in_orig = vol_ct(:,:,:,1);
%     vol_ex_orig = vol_ct(:,:,:,5);
%     
%     init_size = size(vol_in_orig);
%     vol_in = single(crop_data(vol_in_orig, crop_v));
%     vol_ex = single(crop_data(vol_ex_orig, crop_v));
%     avg_img = single(crop_data(avg_img_orig, crop_v));
%     pts_struct = DIR_get_landmarks_for_the_case(num, landmarks_path);
%     pts_in = pts_struct.extreme.in; 
%     pts_ex = pts_struct.extreme.ex;
%     
%     units = info_mask.PixelDimensions;
%     vol_in_nonor = vol_in;
%     vol_ex_nonor = vol_ex;
%     vol_in = img_thr(vol_in, -1024, 0, 1);
%     vol_ex = img_thr(vol_ex, -1024, 0, 1);
%     avg_img = img_thr(avg_img, -1024, 0, 1);
    
    
end