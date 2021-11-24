% follow the https://github.com/dgoodwin208/Registration/wiki
% First you should download the VLFeat toolbox

clear;

addpath(genpath('../../VAMPIRE_FullDatabase_MHA/'));
addpath(genpath('./scr/'));
addpath(genpath('./pTVreg-master/mutils/My/'));
addpath(genpath('./pTVreg-master/ptv'));

% RefVI_type = 'Galligas';
RefVI_type = 'Xenon';
% RefVI_type = 'DTPA-SPECT';

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
            subject_num = 21;
            num_phase = 10;
        end
    end
end
TREs = zeros(subject_num,1);
TREs_or = zeros(subject_num,1);
SpR = zeros(1, subject_num);


save_path = 'H:/CTVI/Results/Continuous_figure/SD02/';
for num = 3:subject_num
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
        info = mha_read_header(CT_Filename);
        vol_ct_cur = mha_read_volume(info);
        vol_ct = cat(4,vol_ct,vol_ct_cur);
    end
    
    for num_idx = 1:num_phase
    vol_in_orig = double(vol_ct(:,:,:,num_idx));
    vol_ex_orig = double(vol_ct(:,:,:,8));
    vol_in_mask = double(vol_mask(:,:,:,num_idx));
    vol_ex_mask = double(vol_mask(:,:,:,8));
    
    pts_struct = DIR_get_landmarks_for_the_case_02(num, landmarks_path);
    pts_in = pts_struct.extreme.in;
    pts_ex = pts_struct.extreme.ex;

    units = info_mask.PixelDimensions;
    
    % isoPTV reg
    use_refinement = 0 ;
    resize = 1;
    fast_lcc =  1;
    spc_orig = units;
    [pt_errs_phys_or, TRE_phys_or, TREstd_phys_or] = DIR_or(pts_in, pts_ex, spc_orig, []);
    TREs_or(num) = TRE_phys_or;
    
    if resize
        bszv = size(vol_in_orig);
        spc_tmp = [1, 1, 1];
        vol_ex = volresize(vol_ex_orig, round(bszv .* units .* spc_tmp), 1);
        vol_in = volresize(vol_in_orig, round(bszv .* units .* spc_tmp), 1);
        
        avg_mask = imresize3(vol_in_mask, round(bszv .* units .* spc_tmp), 'nearest');
        spc = [1,1,1] ./ spc_tmp;
    end
    d = [10, 10, 5];
    vol_rsz = size(vol_ex);
    crop_v = crop_mask(avg_mask,size(avg_mask),d);
    sup_real_sz = [crop_v(1,2)-crop_v(1,1)+1,crop_v(2,2)-crop_v(2,1)+1,crop_v(3,2)-crop_v(3,1)+1];
    vol_ex = crop_data(vol_ex,crop_v);
    vol_in = crop_data(vol_in,crop_v);
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
    
    savename = [save_path,'subject_',num2str(num),num2str(num_idx),'.mat'];
    save(savename,'Tptv');
    [~,Tptv] = uncrop_data(vol_ex,Tptv,crop_v,vol_rsz);
    if resize
        Tptv_rsz = cat(4, volresize(Tptv(:,:,:,1), bszv), volresize(Tptv(:,:,:,2), bszv), volresize(Tptv(:,:,:,3), bszv));
    else
        Tptv_rsz = Tptv;
    end
    Tptv_rsz_all = Tptv_rsz;
    % move points and measure TRE
    [pt_errs_phys, pts_moved_pix, TRE_phys, TREstd_phys] = DIR_movepoints(pts_in, pts_ex, Tptv_rsz_all, spc_orig, []);
    TREs(num) = mean(TRE_phys);
    end
end