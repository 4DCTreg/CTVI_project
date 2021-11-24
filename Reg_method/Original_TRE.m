% This is the tese of the POPI dataset
clear;

addpath(genpath('../pTVreg-master/mutils/My/'));
addpath(genpath('../pTVreg-master/ptv'));
addpath(genpath('../HOMRF_groupwise/MIND-SSC/'));
addpath(genpath('../HOMRF_groupwise/scr/'));
addpath(genpath('../DIRMatlabPackage/'));


isoPTVpath = 'H:/groupwise_registration/POPI_iso_continus/';
LMpath = 'H:/POPI_model/example/4DLandmarks/';
IMpath = 'H:/POPI_model/example/4DCT_MetaImage/';
infpath = [IMpath,'00-P.mhd'];
info = mha_read_header(infpath);
spc_or = info.PixelDimensions;
pts_all_phase = zeros(40,3,10);
koef = repmat(spc_or,40,1);
all_size = info.Dimensions; 
volmov = zeros(all_size);
dvfpath = 'H:/POPI_model/example/4DVF/Parametric/';
obmethod = 'isoPTV';
init_size = info.Dimensions;
ir = init_size(1);
ic = init_size(2);
is = init_size(3);
[Xgrid,Ygrid,Zgrid] = ndgrid(1:ir,1:ic,1:is);
Xgrid = Xgrid.*spc_or(1);
Ygrid = Ygrid.*spc_or(2);
Zgrid = Zgrid.*spc_or(3);
TRE_or = zeros(1,10);
STD_or = zeros(1,10);
filepath = 'H:/POPI_model/example/4DCT_MetaImage/';
    maskpath = 'H:/POPI_model/example/4DMask_metalmage/';
    LMpath = 'H:/POPI_model/example/4DLandmarks/';
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


        init_size = info.Dimensions;
        
        volmov = readrawPOPImeta(filename,init_size);
        volfix = readrawPOPImeta('H:/POPI_model/example/4DCT_MetaImage/10_P.raw',init_size);


        
        
        pts_fix_all = read_POPI_points_file('H:/POPI_model/example/4DLandmarks/case10.txt');
        pts_mov_all = read_POPI_points_file(LMname);
        pts_fix_or = pts_fix_all(1:40,:);
        pts_mov_or = pts_mov_all(1:40,:);

       Knots_n_kalman_rsz_all = zeros(init_size(1),init_size(2),init_size(3),3);

        griddedInterpolantX = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_kalman_rsz_all(:,:,:,1));
        griddedInterpolantY = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_kalman_rsz_all(:,:,:,2));
        griddedInterpolantZ = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_kalman_rsz_all(:,:,:,3));
        displ(:,1) = griddedInterpolantX(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displ(:,2) = griddedInterpolantY(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displ(:,3) = griddedInterpolantZ(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        warpedLMs = displ+pts_fix_or;

        tre = mean(sqrt(sum(((warpedLMs-pts_mov_or)).^2,2)));
        TRE_or(phase) = tre;
        STD_or(phase) = std(sqrt( sum(( (warpedLMs-pts_mov_or)).^2, 2) ));
      
    end


