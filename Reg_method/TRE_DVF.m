% This code is test parametirc and non-parametric method'
clear;

addpath(genpath('../pTVreg-master/mutils/My/'));
addpath(genpath('../pTVreg-master/ptv'));
addpath(genpath('../HOMRF_groupwise/MIND-SSC/'));
addpath(genpath('../HOMRF_groupwise/scr/'));

dvfpath = 'H:/POPI_model/example/4DVF/Parametric/';
T_ob = [];
for phase = 1:10
    if phase ==10
        VFname = [dvfpath,'VF10.vf'];
    else
        VFname = [dvfpath,'VF1',num2str(phase),'.vf'];
    end
    T_es_cur = readvf(VFname);
    T_ob(:,:,:,:,phase) = T_es_cur.*2;

end
init_size = size(T_ob(:,:,:,1,1));

resize = 1;
isoPTVpath = 'H:/groupwise_registration/POPI_iso_continus/';
LMpath = 'H:/POPI_model/example/4DLandmarks/';
IMpath = 'H:/POPI_model/example/4DCT_MetaImage/';
infpath = [IMpath,'00-P.mhd'];
info = mha_read_header(infpath);
spc_or = info.PixelDimensions;
pts_all_phase = zeros(40,3,10);
koef = repmat(spc_or,40,1);
all_size = info.Dimensions;



for phase = 1:10
    if phase == 10
        LMname = [LMpath,'case00.txt'];
    else
        LMname = [LMpath,'case',num2str(phase),'0.txt'];
    end
    cur_pts = read_POPI_points_file(LMname);
    pts_all_phase(:,:,phase) = round(cur_pts(1:40,:)./koef);
    
end


 
    

    TRE_homrf = zeros(10,1);
STD = zeros(10,1);

    
    for phase = 1:10
        spc_orig = spc_or;
         if phase == 10
            LMname = [LMpath,'case00.txt'];
         else
            LMname = [LMpath,'case',num2str(phase),'0.txt'];
        end

        



        pts_fix_all = read_POPI_points_file('H:/POPI_model/example/4DLandmarks/case10.txt');
        pts_mov_all = read_POPI_points_file(LMname);
        pts_fix_or = pts_fix_all(1:40,:);
        pts_mov_or = pts_mov_all(1:40,:);

       
        ir = init_size(1);
        ic = init_size(2);
        is = init_size(3);
        nx = [1:ir]*2-1;
        ny = [1:ic]*2-1;
        nz = [1:is]*2-1;
        [NX, NY, NZ] = ndgrid(nx,ny,nz);

        
        

        T_ob_cur = T_ob(:,:,:,:,phase);
        HOMRFgriddedInterpolantX = griddedInterpolant(NX,NY,NZ,T_ob_cur(:,:,:,1));
        HOMRFgriddedInterpolantY = griddedInterpolant(NX,NY,NZ,T_ob_cur(:,:,:,2));
        HOMRFgriddedInterpolantZ = griddedInterpolant(NX,NY,NZ,T_ob_cur(:,:,:,3));
        displhomrf(:,1) = HOMRFgriddedInterpolantX(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displhomrf(:,2) = HOMRFgriddedInterpolantY(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displhomrf(:,3) = HOMRFgriddedInterpolantZ(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        warpedLMsHOMRF = displhomrf+pts_fix_or;

        tre = mean(sqrt(sum(((warpedLMsHOMRF-pts_mov_or)).^2,2)));
        std_cur = std(sqrt( sum((  (warpedLMsHOMRF-pts_mov_or)  ).^2, 2) ));
        TRE_homrf(phase) = tre;
        STD(phase) = std_cur;


        
    
        
    end


