% follow the https://github.com/dgoodwin208/Registration/wiki
% First you should download the VLFeat toolbox 

clear;

addpath(genpath('../../VAMPIRE_FullDatabase_MHA/'));
addpath(genpath('./scr/'));
addpath(genpath('./pTVreg-master/mutils/My/'));
addpath(genpath('./pTVreg-master/ptv'));

RefVI_type = 'Galligas';
% RefVI_type = 'DTPA-SPECT';
% RefVI_type = 'DTPA-SPECT';

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
TREs = zeros(25,1);
TREs_or = zeros(25,1);
SpR = zeros(1, 25);
landmarks_SIFT = cell(1,25);

save_path = 'H:/CTVI/Results/SD_01_in_to_ex_iso_111/';
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
    GT_mask_filename = [Ref_path,'VentMask.mha'];
    info_mask_GT = mha_read_header(GT_mask_filename);
    GT_mask = mha_read_volume(info_mask_GT);
    
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
    crop_v_GT = crop_mask(GT_mask,info_mask_GT,d);
    GT_mask_crop = single(crop_data(GT_mask, crop_v_GT));
    GT_img_crop = single(crop_data(GT_img, crop_v_GT));
    vol_in_orig = vol_ct(:,:,:,1);
    vol_ex_orig = vol_ct(:,:,:,5);
    
    init_size = size(vol_in_orig);
    vol_in = single(crop_data(vol_in_orig, crop_v));
    vol_ex = single(crop_data(vol_ex_orig, crop_v));
    avg_img = single(crop_data(avg_img_orig, crop_v));
    
%     show_overlay(avg_img_orig, Glandmarks_SIFTT_img, GT_mask, 30, 'autumn') 
    vol_in_nonor = vol_in;
    vol_ex_nonor = vol_ex;
    vol_in = img_thr(vol_in, -1024, 0, 1);
    vol_ex = img_thr(vol_ex, -1024, 0, 1);
    avg_img = img_thr(avg_img, -1024, 0, 1);


    units = info_mask.PixelDimensions;

    
    % isoPTV reg
    use_refinement = 0 ;
    resize = 1;
    fast_lcc =  1;
    spc_orig = units;

    if resize
        bszv = size(vol_in);
        spc_tmp = [1, 1, 1];
        vol_ex = volresize(vol_ex, round(bszv .* units .* spc_tmp), 1);
        vol_in = volresize(vol_in, round(bszv .* units .* spc_tmp), 1);
        vol_ex_nonor = volresize(vol_ex_nonor, round(bszv .* units .* spc_tmp), 1);
        vol_in_nonor = volresize(vol_in_nonor, round(bszv .* units .* spc_tmp), 1);
        avg_img = volresize(avg_img, round(bszv .* units .* spc_tmp), 1);
        spc = [1,1,1] ./ spc_tmp;
    end
    
    % configure registration
    opts = [];
    opts.loc_cc_approximate = fast_lcc;
    if use_refinement
        if resize
            opts.grid_spacing = [4, 4, 4]*2; 
        else
            opts.grid_spacing = [4, 4, 3]*2;  % grid spacing in pixels
        end
        opts.cp_refinements = 1;
    else
        if resize
            opts.grid_spacing = [4, 4, 4]; 
        else
            opts.grid_spacing = [4, 4, 3];  % grid spacing in pixels
        end
        opts.cp_refinements = 0;
    end
    opts.display = 'off';
    opts.k_down = 0.7;
    opts.interp_type = 0;
    opts.metric = 'loc_cc_fftn_gpu';
    opts.metric_param = [1,1,1] * 2.1;

    opts.scale_metric_param = true;
    opts.isoTV = 0.11;
    opts.csqrt = 5e-3;
    opts.spline_order = 1;
    opts.border_mask = 5;
    opts.max_iters =  80;
    opts.check_gradients = 100*0;
    opts.pix_resolution = spc;
    [voldef, Tptv, Kptv] = ptv_register(vol_in, vol_ex, opts);

    savename = [save_path,'subject_',num2str(num),'.mat'];
    save(savename,'Tptv');  
end