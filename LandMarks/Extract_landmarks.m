% follow the https://github.com/dgoodwin208/Registration/wiki
% First you should download the VLFeat toolbox 

clear;

% addpath(genpath('E:\pape_code\prepocessing'));
% addpath(genpath('./3D_SIFT_pkg/'))
addpath(genpath('../../VAMPIRE_FullDatabase_MHA/'));
addpath(genpath('../scr/'));
addpath(genpath('G:/SIFT3D/SIFT3D 1.4.5/'));

RefVI_type = 'Galligas';
% RefVI_type = 'DTPA-SPECT';
% RefVI_type = 'DTPA-SPECT';

if strcmp(RefVI_type,'Galligas')
    file_path = 'H:/VAMPIRE_FullDatabase_MHA/Study01_Galligas-PET/';
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
    vent_Ref = mha_read_volume(info);
    vol_ct = [];
    if strcmp(RefVI_type,'Xenon')
        avg_img = [];
    else
        avg_img_name = [CT_path,'AverageImage.mha'];
        info = mha_read_header(avg_img_name);
        avg_img = mha_read_volume(info);
        avg_mask_name = [CT_path,'AverageMask.mha'];
        info_mask = mha_read_header(avg_mask_name);
        avg_mask = mha_read_volume(info_mask);
    end
    coronal = [];
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
    sz_ct = info.Dimensions;
    d = [10, 10, 5];
    crop_v = crop_mask(avg_mask,info_mask,d);
    vol_in = vol_ct(:,:,:,1);
    vol_ex = vol_ct(:,:,:,4);
    vol_in = single(crop_data(vol_in, crop_v));
    vol_ex = single(crop_data(vol_ex, crop_v));

    units = info_mask.PixelDimensions;
    % Detect keypoints
    keys_in = detectSift3D(vol_in, 'units', units);
    keys_ex = detectSift3D(vol_ex, 'units', units);
    % Extract descriptors
    [desc_in, coords_in] = extractSift3D(keys_in);
    [desc_ex, coords_ex] = extractSift3D(keys_ex);
    correspondences = vl_ubcmatch(desc_in',desc_ex');
  
end

