% follow the https://github.com/dgoodwin208/Registration/wiki
% First you should download the VLFeat toolbox

clear;
addpath(genpath('../CTVI_data/VAMPIRE_FullDatabase_MHA/'));
addpath(genpath('./scr/'));
addpath(genpath('./pTVreg-master/mutils/My/'));
addpath(genpath('./pTVreg-master/ptv'));

% RefVI_type = 'Galligas';
% RefVI_type = 'Xenon';
RefVI_type = 'DTPA-SPECT';

if strcmp(RefVI_type,'Galligas')
    file_path = 'D:/CTVI_data/VAMPIRE_FullDatabase_MHA/Study01_Galligas-PET/';
    subject_num = 25;
    num_phase = 5;
    spc_CT_all = zeros(subject_num,3);
    size_CT_all = zeros(subject_num,3);
    spc_GT_all = zeros(subject_num,3);
    size_GT_all = zeros(subject_num,3);
    max_value_all = zeros(subject_num);
else
    if strcmp(RefVI_type,'Xenon')
        file_path = 'D:/CTVI_data/VAMPIRE_FullDatabase_MHA/Study02_Xenon-CT/';
        subject_num = 4;
        num_phase = 8;
    else
        if strcmp(RefVI_type,'DTPA-SPECT')
            file_path = 'D:/CTVI_data/VAMPIRE_FullDatabase_MHA/Study03_DTPA-SPECT/';
            subject_num = 21;
            num_phase = 10;
            spc_CT_all = zeros(subject_num,3);
            size_CT_all = zeros(subject_num,3);
            spc_GT_all = zeros(subject_num,3);
            size_GT_all = zeros(subject_num,3);
        end
    end
end

for num = 1:subject_num
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
    max_value_all(num) = max(GT_img(:));
    GT_mask_filename = [Ref_path,'VentMask.mha'];
    info_mask_GT = mha_read_header(GT_mask_filename);
    GT_mask = mha_read_volume(info_mask_GT);
    
    spc_GT_all(num, :) = info.PixelDimensions;
    size_GT_all(num,:) = info.Dimensions;
    
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
    
    spc_CT_all(num,:) = info.PixelDimensions;
    size_CT_all(num,:) = info.Dimensions;
    
    
end
disp(mean(spc_GT_all))
disp(mean(size_GT_all))
disp(mean(spc_CT_all))
disp(mean(size_CT_all))
