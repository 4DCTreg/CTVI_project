% This is a test of the POPI dataset and based on 'kalman_test.m', where
% HOMRF method using multi-level processing strategy
clear;

addpath(genpath('../isoPTVpackage/'));
addpath(genpath('./MIND-SSC/'));
addpath(genpath('./scr/'));

IMpath = 'H:/POPI_model/example/4DCT_MetaImage/';
maskpath = 'H:/POPI_model/example/4DMask_metalmage/';
LMpath = 'H:/POPI_model/example/4DLandmarks/';
filepath = 'H:/POPI_model/example/4DCT_MetaImage/';
infpath = [IMpath,'00-P.mhd'];
info = mha_read_header(infpath);
spc_or = info.PixelDimensions;
pts_all_phase = zeros(40,3,10);
koef = repmat(spc_or,40,1);
all_size = info.Dimensions; 
volmov = zeros(all_size);

init_size = info.Dimensions;
ir = init_size(1);
ic = init_size(2);
is = init_size(3);
[Xgrid,Ygrid,Zgrid] = ndgrid(1:ir,1:ic,1:is);
Xgrid = Xgrid.*spc_or(1);
Ygrid = Ygrid.*spc_or(2);
Zgrid = Zgrid.*spc_or(3);

for phase = 1:10
    if phase == 10
        LMname = [LMpath,'case00.txt'];
    else
        LMname = [LMpath,'case',num2str(phase),'0.txt'];
    end
    cur_pts = read_POPI_points_file(LMname);
    pts_all_phase(:,:,phase) = round(cur_pts(1:40,:)./koef);
end

min1 = min(min(squeeze(pts_all_phase(:,1,:))));
max1 = max(max(squeeze(pts_all_phase(:,1,:))));
min2 = min(min(squeeze(pts_all_phase(:,2,:))));
max2 = max(max(squeeze(pts_all_phase(:,2,:))));
min3 = min(min(squeeze(pts_all_phase(:,3,:))));
max3 = max(max(squeeze(pts_all_phase(:,3,:))));
d = [10, 10, 5];
min_max = [min1,max1;min2,max2;min3,max3];
crop_v = [ max(1, min_max(1,1) - d(1)), min(size(volmov, 1), min_max(1,2) + d(1)); ... 
               max(1, min_max(2,1) - d(2)), min(size(volmov, 2), min_max(2,2) + d(2)); ... 
               max(1, min_max(3,1) - d(3)), min(size(volmov, 3), min_max(3,2) + d(3));];
           


    TRE_homrf = zeros(10,1);
    STD_homrf = zeros(10,1);


    for phase = 1:10
        
        if phase == 10
            Movmaskname = [maskpath,'00-air-body-lungs.raw'];
            filename = [filepath,'00_P.raw'];
            LMname = [LMpath,'case00.txt'];
        else
            Movmaskname = [maskpath,num2str(phase),'0-air-body-lungs.raw'];
            filename = [filepath,num2str(phase),'0_P.raw'];
            LMname = [LMpath,'case',num2str(phase),'0.txt'];
        end
        Fixmaskname = [maskpath,'10-air-body-lungs.raw'];
        mov_mask = readrawPOPImask(Movmaskname,init_size);
        fix_mask = readrawPOPImask(Fixmaskname,init_size);
        mov_mask(mov_mask~=2) = 0;
        fix_mask(fix_mask~=2) = 0;
        mov_mask(mov_mask==2) = 1;
        fix_mask(fix_mask==2) = 1;
        spc = spc_or;
        init_size = info.Dimensions;
        
        volmov = readrawPOPImeta(filename,init_size);
        volfix = readrawPOPImeta('H:/POPI_model/example/4DCT_MetaImage/10_P.raw',init_size);
        volmov = volmov+1024;
        volfix = volfix+1024;
        volmov = single(img_thr(volmov, min(volmov(:)), max(volmov(:)), 1));
        volfix = single(img_thr(volfix, min(volfix(:)), max(volfix(:)), 1));

         
        pts_fix_all = read_POPI_points_file('H:/POPI_model/example/4DLandmarks/case10.txt');
        pts_mov_all = read_POPI_points_file(LMname);
        pts_fix_or = pts_fix_all(1:40,:);
        pts_mov_or = pts_mov_all(1:40,:);
        pts_fix = round(pts_fix_or./koef);
        pts_mov = round(pts_mov_or./koef);

        volmov = crop_data(volmov, crop_v);
        volfix = crop_data(volfix, crop_v);
        mov_mask = crop_data(mov_mask,crop_v);
        fix_mask = crop_data(fix_mask,crop_v);

        useTop = 0;
        resize =1;
        if resize==1
            bszv = size(volmov);
            spc_tmp = [1, 1, 1];
            volfix = volresize(volfix, round(bszv .* spc .* spc_tmp), 1);
            volmov = volresize(volmov, round(bszv .* spc .* spc_tmp), 1);
            mov_mask = imresize3(mov_mask, round(bszv .* spc .* spc_tmp), 'nearest');
            fix_mask = imresize3(fix_mask, round(bszv .* spc .* spc_tmp), 'nearest');
            spc = [1,1,1] ./ spc_tmp;
        end
        parameter = homrf_get_POPI_parameter_multi_level(useTop,spc_or);
        
        cur_grid_space = parameter.grid_space(5,:);
        
        [O_trans_X,O_trans_Y,O_trans_Z,dx,dy,dz]=make_control_grid_3D_Odd(cur_grid_space,size(volmov));
              
        griddim = [size(O_trans_X,1),size(O_trans_X,2),size(O_trans_X,3)];
        
        if phase == 1
            Knots_n = zeros(size(volmov,1),size(volmov,2),size(volmov,3),3);
        else
            Knots_n = Knots_n_tot;
        end
 
        [Knots_n_tot,Knots_n] = register_homrf_kalman_POPI_mask_multi_level(parameter,volmov,volfix,mov_mask,fix_mask,...
            init_size,pts_fix_or,pts_mov_or,Knots_n,bszv,crop_v); 
    
        Knots_n_pre = Knots_n;

        if parameter.resize ==1
            Tptv_rsz = cat(4, volresize(Knots_n_tot(:,:,:,1), bszv), ...
       volresize(Knots_n_tot(:,:,:,2), bszv), volresize(Knots_n_tot(:,:,:,3), bszv));  
            voldef = volresize(volmov,bszv);
        end
        [~, Knots_n_rsz] = uncrop_data(voldef, Tptv_rsz, crop_v, init_size);

        HOMRFgriddedInterpolantX = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_rsz(:,:,:,1));
        HOMRFgriddedInterpolantY = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_rsz(:,:,:,2));
        HOMRFgriddedInterpolantZ = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_rsz(:,:,:,3));
        displhomrf(:,1) = HOMRFgriddedInterpolantX(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displhomrf(:,2) = HOMRFgriddedInterpolantY(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displhomrf(:,3) = HOMRFgriddedInterpolantZ(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        warpedLMsHOMRF = displhomrf+pts_fix_or;

        tre = mean(sqrt(sum(((warpedLMsHOMRF-pts_mov_or)).^2,2)));
        TRE_homrf(phase) = tre;
        STD_homrf(phase) = std(sqrt( sum(( (warpedLMsHOMRF-pts_mov_or)).^2, 2) ));
   
    end
