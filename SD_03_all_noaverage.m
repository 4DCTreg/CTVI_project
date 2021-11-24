% follow the https://github.com/dgoodwin208/Registration/wiki
% First you should download the VLFeat toolbox

clear;

addpath(genpath('../../VAMPIRE_FullDatabase_MHA/'));
addpath(genpath('./scr/'));
addpath(genpath('./pTVreg-master/mutils/My/'));
addpath(genpath('./pTVreg-master/ptv'));

% RefVI_type = 'Galligas';
% RefVI_type = 'Xenon';
RefVI_type = 'DTPA-SPECT';
% VI_metric = 'DIR_Jac';
% VI_metric = 'DIR_HU';
VI_metric = 'Supervoxel';
% VI_metric = 'Supervoxel';
% VI_metric = 'Average';

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
            CTVI_path = 'H:/CTVI/Results/CTVI_results/SD_03/';
            subject_num = 21;
            num_phase = 10;
        end
    end
end
TREs = zeros(1, subject_num);
TREs_or = zeros(1, subject_num);
SpR = zeros(1, subject_num);
Dice_all = zeros(subject_num,5);


in_to_ex_path = 'H:/CTVI/Results/SD_03_in_to_ex_iso_multi_sup/';
ex_to_avg_path = 'H:/CTVI/Results/SD_03_ex_to_avg_sup_multi/';
for num = 1:subject_num
    CTVI_name = [CTVI_path,'iso',VI_metric,'_Subject_0',num2str(num),'.mat'];
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
    
    pts_struct = DIR_get_landmarks_for_the_case_03(num, landmarks_path);
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
        avg_mask = imresize3(avg_mask_orig, round(bszv .* units .* spc_tmp), 'nearest');
        spc = [1,1,1] ./ spc_tmp;
    end
    d = [10, 10, 5];
    vol_rsz = size(vol_ex);
    crop_v = crop_mask(avg_mask,size(avg_mask),d);
    sup_real_sz = [crop_v(1,2)-crop_v(1,1)+1,crop_v(2,2)-crop_v(2,1)+1,crop_v(3,2)-crop_v(3,1)+1];
    vol_ex_ROI = crop_data(vol_ex,crop_v);
    vol_in_ROI = crop_data(vol_in,crop_v);
    %     vol_ex_ROI = vol_ex;
    %     vol_in_ROI = vol_in;
    n_seeds = 5000;
    distance = round(nthroot(prod(sup_real_sz)/n_seeds,3));
    h_distance = floor(distance/2);
    [new_crop,min_h_distance] = mulit_sup_corp(crop_v,vol_rsz,h_distance);
    % offset_range
    [offset_Y,offset_X,offset_Z]= meshgrid(-min_h_distance(2):min_h_distance(2),-min_h_distance(1):min_h_distance(1),...
        -min_h_distance(3):min_h_distance(3));
    
    
    %     show_overlay(avg_img_orig, Glandmarks_SIFTT_img, GT_mask, 30, 'autumn')
    %     vol_in_nonor = vol_in;
    %     vol_ex_nonor = vol_ex;
    %     vol_in = img_thr(vol_in, -1024, 0, 1);
    %     vol_ex = img_thr(vol_ex, -1024, 0, 1);
    
    
    in_to_ex_name = [in_to_ex_path,'subject_',num2str(num),'.mat'];
    load(in_to_ex_name);
    
    if strcmp(VI_metric,'Supervoxel')
                [vent_jac_interp] = Ventilation_Jac_coef(vol_ex_ROI, Tptv, spc);
        
        vent_img_multi = zeros(sup_real_sz(1),sup_real_sz(2),sup_real_sz(3));
        [vent_HU] = Ventilation_HU_SD_03(vol_ex_ROI,vol_in_ROI, Tptv, spc);
        for multi_layer = 1:1
            off_cur_index = unidrnd(numel(offset_Y));
            off_X_cur = offset_X(off_cur_index);
            off_Y_cur = offset_Y(off_cur_index);
            off_Z_cur = offset_Z(off_cur_index);
            sup_crop_v_cur = zeros(3,2);
            if off_X_cur<0
                sup_crop_v_cur(1,1) = crop_v(1,1)+off_X_cur;
                sup_crop_v_cur(1,2) = crop_v(1,2);
            else
                sup_crop_v_cur(1,2) = crop_v(1,2)+off_X_cur;
                sup_crop_v_cur(1,1) = crop_v(1,1);
            end
            if off_Y_cur<0
                sup_crop_v_cur(2,1) = crop_v(2,1)+off_Y_cur;
                sup_crop_v_cur(2,2) = crop_v(2,2);
            else
                sup_crop_v_cur(2,2) = crop_v(2,2)+off_Y_cur;
                sup_crop_v_cur(2,1) = crop_v(2,1);
            end
            
            if off_Z_cur<0
                sup_crop_v_cur(3,1) = crop_v(3,1)+off_Z_cur;
                sup_crop_v_cur(3,2) = crop_v(3,2);
            else
                sup_crop_v_cur(3,2) = crop_v(3,2)+off_Z_cur;
                sup_crop_v_cur(3,1) = crop_v(3,1);
            end
            
            %         vol_in_sup_cur = crop_data(vol_in,sup_crop_v_cur);
            vol_ex_sup_cur = crop_data(vol_ex,sup_crop_v_cur);
            multi_sup_sz = size(vol_ex_sup_cur);
            sup_crop_v_cur_recover = zeros(3,2);
            if off_X_cur<0
                sup_crop_v_cur_recover(1,1) = -off_X_cur+1;
                sup_crop_v_cur_recover(1,2) = multi_sup_sz(1);
            else
                sup_crop_v_cur_recover(1,2) = multi_sup_sz(1)-off_X_cur;
                sup_crop_v_cur_recover(1,1) = 1;
            end
            if off_Y_cur<0
                sup_crop_v_cur_recover(2,1) = -off_Y_cur+1;
                sup_crop_v_cur_recover(2,2) = multi_sup_sz(2);
            else
                sup_crop_v_cur_recover(2,2) = multi_sup_sz(2)-off_Y_cur;
                sup_crop_v_cur_recover(2,1) = 1;
            end
            
            if off_Z_cur<0
                sup_crop_v_cur_recover(3,1) = -off_Z_cur+1;
                sup_crop_v_cur_recover(3,2) = multi_sup_sz(3);
            else
                sup_crop_v_cur_recover(3,2) = multi_sup_sz(3)-off_Z_cur;
                sup_crop_v_cur_recover(3,1) = 1;
            end
            [SLIC_Labels_3D,~] = superpixels3(vol_ex_sup_cur,n_seeds);
            SLIC_Labels_3D = crop_data(SLIC_Labels_3D,sup_crop_v_cur_recover);
            labels_index = unique(SLIC_Labels_3D);
            Numlabels = size(labels_index,1);
            
            vol_in_nonor = vol_in_ROI;
            vol_ex_nonor = vol_ex_ROI;
            
            index = 900;
            %             [V_before,C_p_ex,supvoxel_sz] = Supervoxel_Volume(SLIC_Labels_3D,labels_index);
            
            % CTVI img at exhale phase
           [V_before,V_after,N_in_after,C_p_in] = Supervoxel_in(SLIC_Labels_3D,Tptv,spc,labels_index);
            %             vent_img_cur = Ventilation_Supervoxel_HU(vol_ex_nonor,vol_in_nonor,SLIC_Labels_3D,...
            %     labels_index,Tptv,spc);
            [vent_img_cur] = Ventilation_Supervoxel_SD_03(vol_ex_nonor,vol_in_nonor,SLIC_Labels_3D,vent_HU,V_after,...
                C_p_in,labels_index,Tptv,spc);
%             [vent_img_cur] = Ventilation_Supervoxel_no_linear_19(vol_ex_nonor,vol_in_nonor,vent_jac_interp,SLIC_Labels_3D,...
%                         V_after, C_p_in,labels_index,Tptv,spc);
            vent_img_multi = vent_img_cur+vent_img_multi;
        end
        vent_img = vent_img_multi;
    else
        %         [V_after,C_p_in] = Supervoxel_in_forward_backward(SLIC_Labels_3D,Tptv,spc,labels_index);
        if strcmp(VI_metric,'DIR_Jac')
            [vent_img] = Ventilation_Jac(vol_ex_ROI, Tptv, spc);
        else
            if strcmp(VI_metric,'DIR_HU')
                [vent_img] = Ventilation_HU(vol_ex_ROI,vol_in_ROI, Tptv, spc);
            else
                if strcmp(VI_metric,'Average')
                    vent_img = Ventilation_HU_average(vol_ct);
                end
            end
        end
        
    end
%     save(CTVI_name,'vent_img');
    if strcmp(VI_metric,'Average')
        vent_img_mask = vent_img.*double(avg_mask_orig);
        
        GT_img_sz = size(GT_img);
        
        avg_img_HU_rsz = volresize(avg_img_orig,GT_img_sz, 1);
        avg_img_HU_rsz(find(avg_img_HU_rsz==0)) = 0.0001;
        
        % resize CTVI
        vent_img_mask = volresize(vent_img_mask, GT_img_sz, 1);
        vent_img_mask = vent_img_mask+0.01;
        
        GT_img = GT_img + 0.01;
        
        
        % median filter
                vent_img_avg_filter = medfilt3(vent_img_mask,[3,3,3]);
                GT_img_filter = medfilt3(GT_img,[3,3,3]);

        index_slice = 256;
        figure
        show_overlay_sagittal_SD03(avg_img_orig, GT_img_filter, info_GT, GT_mask, index_slice,'hot')
        figure
        show_overlay_sagittal_CTVI_SD_03(avg_img_orig, vent_img_avg_filter, info_GT, GT_mask, index_slice,'hot')
        
        vent_img_avg_filter = vent_img_mask;
        GT_img_filter = GT_img;
        %     GT_img_filter = GT_img;
        
        vent_jac = double(vent_img_avg_filter).*double(GT_mask);
        vent_GT = double(GT_img_filter).*double(GT_mask);
        Dice_cur = Dice_coef(vent_jac,vent_GT,GT_img_sz);
        Dice_all(num,:) = Dice_cur;
        %     show_overlay(avg_img_orig, vent_img_avg_filter, GT_mask, 30, 'autumn');
        %     show_overlay(avg_img_orig, GT_img_filter, GT_mask, 30, 'autumn');
        Compare_vent_jac = reshape(vent_jac,[GT_img_sz(1)*GT_img_sz(2)*GT_img_sz(3),1]);
        Compare_vent_GT = reshape(vent_GT,[GT_img_sz(1)*GT_img_sz(2)*GT_img_sz(3),1]);
        Compare_vent_jac(find(Compare_vent_jac==0)) = [];
        Compare_vent_GT(find(Compare_vent_GT==0)) = [];
        SpR(num) = corr(Compare_vent_jac,Compare_vent_GT,'type','Spearman');
    else
        
        
        [~,Tptv] = uncrop_data(vol_ex_ROI,Tptv,crop_v,vol_rsz);
        if resize
            Tptv_rsz = cat(4, volresize(Tptv(:,:,:,1), bszv), volresize(Tptv(:,:,:,2), bszv), volresize(Tptv(:,:,:,3), bszv));
        else
            Tptv_rsz = Tptv;
        end
        Tptv_rsz_all = Tptv_rsz;
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
        [vent_img_avg, Tptv_all] = uncrop_data(vent_img_avg, T_pix_ex_to_avg, crop_v, vol_rsz);
        if resize
            %         Tptv_rsz = cat(4, volresize(T_pix_ex_to_avg(:,:,:,1), bszv), volresize(T_pix_ex_to_avg(:,:,:,2), bszv), volresize(T_pix_ex_to_avg(:,:,:,3), bszv));
            vent_img_avg = volresize(vent_img_avg, bszv);
        end
        
        vent_img_avg = vent_img_avg.*double(avg_mask_orig);
        
        GT_img_sz = size(GT_img);
        
        avg_img_HU_rsz = volresize(avg_img_orig,GT_img_sz, 1);
        avg_img_HU_rsz(find(avg_img_HU_rsz==0)) = 0.0001;
        
        % resize CTVI
        vent_img_avg = volresize(vent_img_avg, GT_img_sz, 1);
        vent_img_avg = vent_img_avg+0.01;
        %     vent_img_avg(find(vent_img_avg==0)) = 0.01;
        % add HU
        %     vent_img_avg = vent_img_avg./(-avg_img_HU_rsz/1000);
        GT_img = GT_img + 0.01;
        
        
        % median filter
            vent_img_avg_filter = medfilt3(vent_img_avg,[3,3,3]);
%         vent_img_avg_filter = vent_img_avg;
        % GT_img_filter = GT_img;
            GT_img_filter = medfilt3(GT_img,[3,3,3]);
%         GT_img_filter = GT_img;
%         index_slice = 256;
%         figure
%         show_overlay_sagittal_SD03(avg_img_orig, GT_img_filter, info_GT, GT_mask, index_slice,'hot')
%         figure
%         show_overlay_sagittal_CTVI_SD_03(avg_img_orig, vent_img_avg_filter, info_GT, GT_mask, index_slice,'hot')
        
        %     vent_img_avg_filter(find(vent_img_avg_filter==0)) = 0.001;
        %     GT_img_filter(find(GT_img_filter ==0)) = 0.001;
        vent_jac = double(vent_img_avg_filter).*double(GT_mask);
        vent_GT = double(GT_img_filter).*double(GT_mask);
        Dice_cur = Dice_coef(vent_jac,vent_GT,GT_img_sz);
        Dice_all(num,:) = Dice_cur;
        %     show_overlay(avg_img_orig, vent_img_avg_filter, GT_mask, 30, 'autumn');
        %     show_overlay(avg_img_orig, GT_img_filter, GT_mask, 30, 'autumn');
        Compare_vent_jac = reshape(vent_jac,[GT_img_sz(1)*GT_img_sz(2)*GT_img_sz(3),1]);
        Compare_vent_GT = reshape(vent_GT,[GT_img_sz(1)*GT_img_sz(2)*GT_img_sz(3),1]);
        Compare_vent_jac(find(Compare_vent_jac==0)) = [];
        Compare_vent_GT(find(Compare_vent_GT==0)) = [];
        SpR(num) = corr(Compare_vent_jac,Compare_vent_GT,'type','Spearman');
    end
end
