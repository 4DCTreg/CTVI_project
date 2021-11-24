% follow the https://github.com/dgoodwin208/Registration/wiki
% First you should download the VLFeat toolbox

clear;

% addpath(genpath('E:\pape_code\prepocessing'));
% addpath(genpath('./3D_SIFT_pkg/'))
addpath(genpath('../../VAMPIRE_FullDatabase_MHA/'));
addpath(genpath('./scr/'));
addpath(genpath('G:/SIFT3D/SIFT3D 1.4.5/'));
addpath(genpath('./MatlabUtilityPack'));
ex_to_avg_path = 'H:/CTVI/Results/SD_01_ex_to_avg/';

addpath(genpath('./pTVreg-master/mutils/My/'));
addpath(genpath('./pTVreg-master/ptv'));

RefVI_type = 'Galligas';
% RefVI_type = 'DTPA-SPECT';
% RefVI_type = 'DTPA-SPECT';
VI_metric = 'DIR_HU';
% VT_metric = 'DIR_Jac';

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
landmarks_SIFT = cell(1,25);
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
    
    vol_in_nonor = vol_in;
    vol_ex_nonor = vol_ex;
    vol_in = img_thr(vol_in, -1024, 0, 1);
    vol_ex = img_thr(vol_ex, -1024, 0, 1);
    avg_img = img_thr(avg_img, -1024, 0, 1);
    
    % load landmarks
    pts_struct = DIR_get_landmarks_for_the_case(num, landmarks_path);
    pts_in = pts_struct.extreme.in;
    pts_ex = pts_struct.extreme.ex;
    
    units = info_mask.PixelDimensions;
    % Automatic landmarks extraction
    %     [pts_in,pts_ex] = SIFT_extract_landmarks(vol_in, vol_ex, units, crop_v, avg_mask);
    %     landmarks_SIFT{num} = cat(3, pts_in, pts_ex);
    
    % isoPTV reg
    use_refinement = 0 ;
    resize = 1;
    fast_lcc =  1;
    spc_orig = units;
    [pt_errs_phys_or, TRE_phys_or, TREstd_phys_or] = DIR_or(pts_in, pts_ex, spc_orig, []);
    TREs_or(num) = TRE_phys_or;
    if resize
        bszv = size(vol_in);
        spc_tmp = [1, 1, 0.5];
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
    
    % CTVI img at exhale phase
    if strcmp(VI_metric,'DIR_Jac')
        [vent_img] = Ventilation_Jac(vol_ex, Tptv, spc_orig);
    else
        if strcmp(VI_metric,'DIR_HU')
            [vent_img] = Ventilation_HU(vol_ex_nonor,vol_in_nonor, Tptv, spc);
        end
    end
    
    if resize
        Tptv_rsz = cat(4, volresize(Tptv(:,:,:,1), bszv), volresize(Tptv(:,:,:,2), bszv), volresize(Tptv(:,:,:,3), bszv));
        voldef = volresize(voldef, bszv);
    else
        Tptv_rsz = Tptv;
    end
    [~, Tptv_rsz_all] = uncrop_data(voldef, Tptv_rsz, crop_v, init_size);
    % move points and measure TRE
    [pt_errs_phys, pts_moved_pix, TRE_phys, TREstd_phys] = DIR_movepoints(pts_in, pts_ex, Tptv_rsz_all, spc_orig, []);
    TREs(num) = mean(TRE_phys);
    
    % match the exhale image to average image
    ex_to_avg_name = [ex_to_avg_path,'subject_',num2str(num),'.mat'];
    load(ex_to_avg_name);
    %     [voldef, Tptv_co, Kptv_co] = ptv_register(vol_ex, avg_img, opts);
    %     T_pix = conv_3d_T_from_phys_to_pix(Tptv_co, spc_orig);
    % get the ventilation image at average phase
    vent_img_avg = ptv_deform(vent_img, T_pix_ex_to_avg, 0);
    if resize
        Tptv_rsz = cat(4, volresize(T_pix_ex_to_avg(:,:,:,1), bszv), volresize(T_pix_ex_to_avg(:,:,:,2), bszv), volresize(T_pix_ex_to_avg(:,:,:,3), bszv));
        vent_img_avg = volresize(vent_img_avg, bszv);
    end
    [vent_img_avg, Tptv_all] = uncrop_data(vent_img_avg, Tptv_rsz, crop_v, init_size);
    vent_img_avg = vent_img_avg.*double(avg_mask);
    
    GT_img_sz = size(GT_img);
    
    avg_img_HU_rsz = volresize(avg_img_orig,GT_img_sz, 1);
    avg_img_HU_rsz(find(avg_img_HU_rsz==0)) = 0.0001;
    
    % resize CTVI
    vent_img_avg = volresize(vent_img_avg, GT_img_sz, 1);
    vent_img_avg = vent_img_avg+0.01;
    %     vent_img_avg(find(vent_img_avg==0)) = 0.01;
    % add HU
    %     vent_img_avg = vent_img_avg./(-avg_img_HU_rsz/1000);
    GT_img = GT_img+0.01;
    %     GT_img(find(GT_img==0)) = 0.01;
    %debug filter
    %     filtername = ['H:\CTVI\Results\2020_6_18_filter_test\','Subject_',num2str(num),'.mat'];
    %     save(filtername, 'vent_img_avg', 'GT_img');
    
    % filter choose
    % mean filter
    %     H_mean = ones(5,5,3)/(5*5*3);
    %     vent_img_avg_filter = imfilter(vent_img_avg,H_mean);
    %     GT_img_filter = imfilter(GT_img,H_mean);
    % gaussian filter
    %     vent_img_avg_filter = imgaussfilt3(vent_img_avg,0.8,'FilterSize',[5,5,3]);
    %     GT_img_filter = imgaussfilt3(GT_img, 0.8, 'FilterSize', [5,5,3]);
    % median filter
    vent_img_avg_filter = medfilt3(vent_img_avg);
    GT_img_filter = GT_img;
    %     vent_img_avg_filter(find(vent_img_avg_filter==0)) = 0.001;
    %     GT_img_filter(find(GT_img_filter ==0)) = 0.001;
    vent_jac = double(vent_img_avg_filter).*double(GT_mask);
    vent_GT = double(GT_img_filter).*double(GT_mask);
    %     show_overlay(avg_img_orig, vent_img_avg_filter, GT_mask, 30, 'autumn');
    %     show_overlay(avg_img_orig, GT_img_filter, GT_mask, 30, 'autumn');
    Compare_vent_jac = reshape(vent_jac,[GT_img_sz(1)*GT_img_sz(2)*GT_img_sz(3),1]);
    Compare_vent_GT = reshape(vent_GT,[GT_img_sz(1)*GT_img_sz(2)*GT_img_sz(3),1]);
    Compare_vent_jac(find(Compare_vent_jac==0)) = [];
    Compare_vent_GT(find(Compare_vent_GT==0)) = [];
    SpR(num) = corr(Compare_vent_jac,Compare_vent_GT,'type','Spearman');
end
[R_tre_sp,p_tre_sp] = corrcoef(SpR,TREs');

save results_iso_dir_HU.mat landmarks_SIFT TREs_or TREs SpR

