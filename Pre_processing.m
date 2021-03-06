% follow the https://github.com/dgoodwin208/Registration/wiki
% First you should download the VLFeat toolbox

clear;
addpath(genpath('../CTVI_data/VAMPIRE_FullDatabase_MHA/'));
addpath(genpath('scr/'));
addpath(genpath('pTVreg-master/mutils/My/'));
addpath(genpath('pTVreg-master/ptv'));
addpath(genpath('Utilize/'));

RefVI_type = 'Galligas';
% RefVI_type = 'DTPA-SPECT';
% RefVI_type = 'DTPA-SPECT';

if strcmp(RefVI_type,'Galligas')
    file_path = 'D:/CTVI_data/VAMPIRE_FullDatabase_MHA/Study01_Galligas-PET/';
    subject_num = 25;
    num_phase = 5;
else
    if strcmp(RefVI_type,'Xenon')
        file_path = 'H:/VAMPIRE_FullDatabase_MHA/Study02_Xenon-CT/';
        subject_num = 4;
        num_phase = 8;
    else
        if strcmp(RefVI_type,'DTPA-SPECT')
            file_path = 'D:/CTVI_data/VAMPIRE_FullDatabase_MHA/Study03_DTPA-SPECT/';
            subject_num = 21;
            num_phase = 10;
        end
    end
end

save_path = 'D:/CTVI_data/Step1_pre_processing/';
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
        avg_mask_orig = mha_read_volume(info_mask);
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
    
    spc_ct = info.PixelDimensions;
    spc_ref = info_mask_GT.PixelDimensions;
    
    
    spc_tmp =  [0.5, 0.5, 0.5];
    
    GT_sz = size(GT_img);
    img_ref = volresize(GT_img, round(GT_sz .* spc_ref .* spc_tmp), 1);
    img_ref_mask =  round(volresize(GT_mask, round(GT_sz .* spc_ref .* spc_tmp), 1));
    bszv = size(avg_mask_orig);
    avg_img = volresize(avg_img_orig, round(bszv .* spc_ct .* spc_tmp), 1);
    avg_mask = round(volresize(avg_mask_orig, round(bszv .* spc_ct .* spc_tmp), 1));
    
    if strcmp(RefVI_type,'Galligas')
        [optimizer,metric] = imregconfig('monomodal');
        tform = imregtform(avg_mask, img_ref_mask, 'rigid', optimizer, metric);
        Rmoving = imref3d(size(avg_mask));
        Rfixed = imref3d(size(img_ref_mask));
        moved_avg_img = imwarp(avg_img,Rmoving,tform,'OutputView',Rfixed, 'SmoothEdges', true);
        moved_avg_mask = imwarp(avg_mask,Rmoving,tform,'OutputView',Rfixed, 'SmoothEdges', true);
    else
        moved_avg_img = avg_img;
        moved_avg_mask = avg_mask;
    end
    vol_cur_rigid_all = [];
    d = [2, 2, 2];
    crop_v = crop_mask(img_ref_mask,size(img_ref_mask),d);
    size_crop_mask = [crop_v(1,2)-crop_v(1,1)+1,crop_v(2,2)-crop_v(2,1)+1, crop_v(3,2)-crop_v(3,1)+1];
    
    % crop for U-net
    r_or = size_crop_mask(1);
    c_or = size_crop_mask(2);
    s_or = size_crop_mask(3);
    
    r_new = ceil(r_or/16)*16;
    c_new = ceil(c_or/16)*16;
    s_new = ceil(s_or/16)*16;
    res_r = r_new - r_or;
    res_c = c_new - c_or;
    res_s = s_new - s_or;
    
    if res_r==0
        r_unet1 = crop_v(1,1);
        r_unet2 = crop_v(1,2);
    else
        r_unet1 = crop_v(1,1) - floor(res_r/2);
        r_unet2 = crop_v(1,2) + ceil(res_r/2);
    end
    
    if res_c==0
        c_unet1 = crop_v(2,1);
        c_unet2 = crop_v(2,2);
    else
        c_unet1 = crop_v(2,1) - floor(res_c/2);
        c_unet2 = crop_v(2,2) + ceil(res_c/2);
    end
    
    if res_s==0
        s_unet1 = crop_v(3,1);
        s_unet2 = crop_v(3,2);
    else
        s_unet1 = crop_v(3,1) - floor(res_s/2)-4;
        if s_unet1<1
            s_unet1= 1;
            s_unet2 = s_unet1 + s_new - 1;
            if s_unet2 > size(img_ref_mask,3)
                s_new_small = floor(s_or/16)*16;
                s_unet2 = s_unet1 + s_new_small - 1;
            end
                
        else
            s_unet2 = crop_v(3,2) + ceil(res_s/2)-4;
        end
    end
    
    crop_unet = [r_unet1,r_unet2;c_unet1,c_unet2;s_unet1,s_unet2];
    for i_idx = 1:5
        CT_cur = vol_ct(:,:,:,i_idx);
        bszv = size(CT_cur);
        spc_tmp = [0.5, 0.5, 0.5];
        vol_cur = volresize(CT_cur, round(bszv .* spc_ct .* spc_tmp), 1);
        if strcmp(RefVI_type,'Galligas')
            vol_cur_rigid = imwarp(vol_cur,Rmoving,tform,'OutputView',Rfixed, 'SmoothEdges', true);
        else
            vol_cur_rigid = vol_cur;
        end
        %         vol_cur_rigid = crop_data(vol_cur_rigid,crop_v);
        vol_cur_rigid = crop_data(vol_cur_rigid,crop_unet);
        vol_cur_rigid = img_thr(vol_cur_rigid, -1024, 100, 1);
        vol_cur_rigid_all = cat(4, vol_cur_rigid_all, vol_cur_rigid);
        spc = [1,1,1] ./ spc_tmp;
    end
    
    
    
    %     gt_img = crop_data(img_ref,crop_v);
    gt_img = crop_data(img_ref,crop_unet);
    
    %     gt_mask = crop_data(img_ref_mask, crop_v);
    gt_mask = crop_data(img_ref_mask, crop_unet);
    
    gt_img = img_thr(gt_img,prctile(gt_img(gt_mask==1),1),prctile(gt_img(gt_mask==1),99),1);
    
    CT_avg_mask = crop_data(round(moved_avg_mask),crop_unet);
    %     CT_avg_mask = crop_data(CT_avg_mask,crop_unet);
    CT_avg_img = crop_data(moved_avg_img, crop_unet);
    %     CT_avg_img = crop_data(CT_avg_img, crop_unet);
    CT_avg_img = img_thr(CT_avg_img, -1024, 100, 1);
    
    if strcmp(RefVI_type,'Galligas')
        gt_img_name = [save_path,'Sub_',num2str(num),'_GT_img','.nii'];
        gt_mask_name = [save_path,'Sub_',num2str(num),'_GT_mask','.nii'];
        CT_avg_img_name = [save_path,'Sub_',num2str(num),'_CT_avg_img','.nii'];
        CT_avg_mask_name = [save_path,'Sub_',num2str(num),'_CT_avg_mask','.nii'];
        CT_1_img = [save_path,'Sub_',num2str(num),'_CT_img_1','.nii'];
        CT_2_img = [save_path,'Sub_',num2str(num),'_CT_img_2','.nii'];
        CT_3_img = [save_path,'Sub_',num2str(num),'_CT_img_3','.nii'];
        CT_4_img = [save_path,'Sub_',num2str(num),'_CT_img_4','.nii'];
        CT_5_img = [save_path,'Sub_',num2str(num),'_CT_img_5','.nii'];
    else
        gt_img_name = [save_path,'Sub_',num2str(num+25),'_GT_img','.nii'];
        gt_mask_name = [save_path,'Sub_',num2str(num+25),'_GT_mask','.nii'];
        CT_avg_img_name = [save_path,'Sub_',num2str(num+25),'_CT_avg_img','.nii'];
        CT_avg_mask_name = [save_path,'Sub_',num2str(num+25),'_CT_avg_mask','.nii'];
        CT_1_img = [save_path,'Sub_',num2str(num+25),'_CT_img_1','.nii'];
        CT_2_img = [save_path,'Sub_',num2str(num+25),'_CT_img_2','.nii'];
        CT_3_img = [save_path,'Sub_',num2str(num+25),'_CT_img_3','.nii'];
        CT_4_img = [save_path,'Sub_',num2str(num+25),'_CT_img_4','.nii'];
        CT_5_img = [save_path,'Sub_',num2str(num+25),'_CT_img_5','.nii'];
    end
    
    gt_img_nii = make_nii(single(gt_img));
    gt_mask_nii = make_nii(single(gt_mask));
    CT_avg_img_nii = make_nii(single(CT_avg_img));
    CT_avg_mask_nii = make_nii(single(CT_avg_mask));
    CT_1_img_nii = make_nii(single(vol_cur_rigid_all(: ,:,:,1)));
    CT_2_img_nii = make_nii(single(vol_cur_rigid_all(: ,:,:,2)));
    CT_3_img_nii = make_nii(single(vol_cur_rigid_all(: ,:,:,3)));
    CT_4_img_nii = make_nii(single(vol_cur_rigid_all(: ,:,:,4)));
    CT_5_img_nii = make_nii(single(vol_cur_rigid_all(: ,:,:,5)));
    
    save_nii(gt_img_nii,gt_img_name);
    save_nii(gt_mask_nii, gt_mask_name);
    save_nii(CT_avg_img_nii,CT_avg_img_name);
    save_nii(CT_avg_mask_nii,CT_avg_mask_name);
    save_nii(CT_1_img_nii,CT_1_img);
    save_nii(CT_2_img_nii,CT_2_img);
    save_nii(CT_3_img_nii,CT_3_img);
    save_nii(CT_4_img_nii,CT_4_img);
    save_nii(CT_5_img_nii,CT_5_img);
    
    
end