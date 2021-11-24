% This is the same version as 'registration_kalman_v1'
clear;

addpath(genpath('../pTVreg-master/mutils/My/'));
addpath(genpath('../pTVreg-master/ptv'));
addpath(genpath('../HOMRF_groupwise/MIND-SSC/'));
addpath(genpath('../HOMRF_groupwise/scr/'));

dvfpath = 'H:/POPI_model/example/4DVF/Parametric/';
T_ob = [];
for phase = 2:9
    
    VFname = [dvfpath,'VF1',num2str(phase),'.vf'];
    T_es_cur = readvf(VFname);
    T_ob(:,:,:,:,phase) = T_es_cur.*2;

end

resize = 1;
isoPTVpath = 'H:/groupwise_registration/POPI_iso_continus/';
LMpath = 'H:/POPI_model/example/4DLandmarks/';
IMpath = 'H:/POPI_model/example/4DCT_MetaImage/';
infpath = [IMpath,'00-P.mhd'];
info = mha_read_header(infpath);
spc_or = info.PixelDimensions;
pts_all_phase = zeros(40,3,9);
koef = repmat(spc_or,40,1);
all_size = info.Dimensions;
volmov = zeros(all_size);
T_es = [];
for phase = 1:9
    LMname = [LMpath,'case',num2str(phase),'0.txt'];
    fieldFilename = [isoPTVpath,'case',num2str(phase),'_0.mat'];
    load(fieldFilename);
    T_es(:,:,:,:,phase) = Tptv;
    cur_pts = read_POPI_points_file(LMname);
    pts_all_phase(:,:,phase) = round(cur_pts(1:40,:)./koef);
    
end

% T_es = [];
% for phase = 1:9
%     
%     nois_cur = rand(235,176,141,3);
%     
%     LMname = [LMpath,'case',num2str(phase),'0.txt'];
%     fieldFilename = [isoPTVpath,'case',num2str(phase),'_0.mat'];
%     load(fieldFilename);
%     T_es(:,:,:,:,phase) = Tptv+nois_cur./3;
%     cur_pts = read_POPI_points_file(LMname);
%     pts_all_phase(:,:,phase) = round(cur_pts(1:40,:)./koef);
%     
% end


    nrsz = size(T_ob);
    nrk = nrsz(1);nck = nrsz(2);nsk  = nrsz(3);
    N=9;
T_ob(:,:,:,:,1) = zeros(nrk,nck,nsk,3,'single');
    Z = zeros(1,N);
    Px = zeros(nrk,nck,nsk,N,'single')+1; 
    Py = zeros(nrk,nck,nsk,N,'single')+1;
    Pz = zeros(nrk,nck,nsk,N,'single')+1;


    Q = 0.25;%W(k)的方差
    R = 0.25;%V(k)的方差

    F=1;
    G=1;
    H=ones(nrk,nck,nsk,'single');
    I=ones(nrk,nck,nsk,'single'); 
    
    xd_kf = zeros(nrk,nck,nsk,N,'single');
    yd_kf = zeros(nrk,nck,nsk,N,'single');
    zd_kf = zeros(nrk,nck,nsk,N,'single');
    xd_kf(:,:,:,1) = single(repmat((1:nrk)',1,nck,nsk));
    yd_kf(:,:,:,1) = single(repmat((1:nck),nrk,1,nsk));
    curk = zeros(nrk,nck,nsk,'single');
    for ly = 1:nsk
        curk(:,:,ly) = ones(nrk,nck)*ly;
    end
    zd_kf(:,:,:,1) = single(curk);
    fix_all = zeros(nrk,nck,nsk,3,'single');
    fix_all(:,:,:,1) = xd_kf(:,:,:,1);
    fix_all(:,:,:,2) = yd_kf(:,:,:,1);
    fix_all(:,:,:,3) = zd_kf(:,:,:,1);
    filepath = 'H:/POPI_model/example/4DCT_MetaImage/';
    TRE_homrf = zeros(10,1);
    TRE_kalman = zeros(10,1);

    
    for phase = 1:9
        spc_orig = spc_or;

        
        init_size = nrsz;
        LMpath = 'H:/POPI_model/example/4DLandmarks/';
        LMname = [LMpath,'case',num2str(phase),'0.txt'];
        pts_fix_all = read_POPI_points_file('H:/POPI_model/example/4DLandmarks/case10.txt');
        pts_mov_all = read_POPI_points_file(LMname);
        pts_fix_or = pts_fix_all(1:40,:);
        pts_mov_or = pts_mov_all(1:40,:);
        bszv = init_size;
       
        ir = init_size(1);
        ic = init_size(2);
        is = init_size(3);
        nx = [1:ir]*2-1;
        ny = [1:ic]*2-1;
        nz = [1:is]*2-1;
        [NX, NY, NZ] = ndgrid(nx,ny,nz);

        
        
%        if resize==1
%             T_es_cur = cat(4, volresize(T_es(:,:,:,1,phase), bszv), volresize(T_es(:,:,:,2,phase), bszv), volresize(T_es(:,:,:,3,phase), bszv));   
%         end
        T_es_cur = T_es(:,:,:,:,phase);
        HOMRFgriddedInterpolantX = griddedInterpolant(NX,NY,NZ,T_es_cur(:,:,:,1));
        HOMRFgriddedInterpolantY = griddedInterpolant(NX,NY,NZ,T_es_cur(:,:,:,2));
        HOMRFgriddedInterpolantZ = griddedInterpolant(NX,NY,NZ,T_es_cur(:,:,:,3));
        displhomrf(:,1) = HOMRFgriddedInterpolantX(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displhomrf(:,2) = HOMRFgriddedInterpolantY(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displhomrf(:,3) = HOMRFgriddedInterpolantZ(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        warpedLMsHOMRF = displhomrf+pts_fix_or;

        tre = mean(sqrt(sum(((warpedLMsHOMRF-pts_mov_or)).^2,2)));
        TRE_homrf(phase) = tre;


        if phase > 1
            es_all = T_es(:,:,:,:,phase);

            es_mov_all = es_all+fix_all;

            move_all = T_ob(:,:,:,:,phase);
            move_all = move_all+fix_all;
            Fx_cur = es_mov_all(:,:,:,1)./xd_kf(:,:,:,phase-1);
            Fy_cur = es_mov_all(:,:,:,2)./yd_kf(:,:,:,phase-1);
            Fz_cur = es_mov_all(:,:,:,3)./zd_kf(:,:,:,phase-1);
            xd_pre = Fx_cur.*xd_kf(:,:,:,phase-1)+Q;
            yd_pre = Fy_cur.*yd_kf(:,:,:,phase-1)+Q;
            zd_pre = Fz_cur.*zd_kf(:,:,:,phase-1)+Q;
            
            p_pre_x = Fx_cur.*Px(:,:,:,phase-1).*Fx_cur+Q;
            p_pre_y = Fy_cur.*Py(:,:,:,phase-1).*Fy_cur+Q;
            p_pre_z = Fz_cur.*Pz(:,:,:,phase-1).*Fz_cur+Q;
            
            Kg_x = p_pre_x.*H.*1./(H.*p_pre_x.*H+R);
            Kg_y = p_pre_y.*H.*1./(H.*p_pre_y.*H+R);
            Kg_z = p_pre_z.*H.*1./(H.*p_pre_z.*H+R);
            
            ex = move_all(:,:,:,1)-H.*xd_pre;
            ey = move_all(:,:,:,2)-H.*yd_pre;
            ez = move_all(:,:,:,3)-H.*zd_pre;
            
            xd_kf(:,:,:,phase) = xd_pre+Kg_x.*ex;
            yd_kf(:,:,:,phase) = yd_pre+Kg_y.*ey;
            zd_kf(:,:,:,phase) = zd_pre+Kg_z.*ez;
            
            Px(:,:,:,phase) = (I-Kg_x.*H).*p_pre_x;
            Py(:,:,:,phase) = (I-Kg_y.*H).*p_pre_y;
            Pz(:,:,:,phase) = (I-Kg_z.*H).*p_pre_z;
        end
        Knots_n_kalman = [];
        
        Knots_n_kalman(:,:,:,1) = xd_kf(:,:,:,phase)-fix_all(:,:,:,1);
        Knots_n_kalman(:,:,:,2) = yd_kf(:,:,:,phase)-fix_all(:,:,:,2);
        Knots_n_kalman(:,:,:,3) = zd_kf(:,:,:,phase)-fix_all(:,:,:,3);
%        if resize==1
%             Knots_n_kalman_rsz = cat(4, volresize(Knots_n_kalman(:,:,:,1), bszv), volresize(Knots_n_kalman(:,:,:,2), bszv), volresize(Knots_n_kalman(:,:,:,3), bszv));   
%         end
            

        Knots_n_kalman_rsz_all = Knots_n_kalman;


    griddedInterpolantX = griddedInterpolant(NX,NY,NZ,Knots_n_kalman_rsz_all(:,:,:,1));
    griddedInterpolantY = griddedInterpolant(NX,NY,NZ,Knots_n_kalman_rsz_all(:,:,:,2));
    griddedInterpolantZ = griddedInterpolant(NX,NY,NZ,Knots_n_kalman_rsz_all(:,:,:,3));
    displ(:,1) = griddedInterpolantX(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
    displ(:,2) = griddedInterpolantY(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
    displ(:,3) = griddedInterpolantZ(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
    warpedLMs = displ+pts_fix_or;

    tre = mean(sqrt(sum(((warpedLMs-pts_mov_or)).^2,2)));
    TRE_kalman(phase) = tre;
    
        
    end


