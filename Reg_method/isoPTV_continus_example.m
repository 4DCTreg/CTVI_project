clear
addpath(genpath('../pTvreg-master/mutils/My/'));
addpath(genpath('../pTvreg-master/ptv'));
use_refinement = 0;
resize = 1;
fast_lcc = 1;
TRE_iso = zeros(10,1);
TRE_me = zeros(10,1);

spc = [0.976562 0.976562 2];
pts_all_phase = zeros(41,3,9);
koef = repmat(spc,41,1);
LMpath = 'H:/POPI_model/example/4DLandmarks/';
IMpath = 'H:/POPI_model/example/4DCT_MetaImage/';
infpath = [IMpath,'00-P.mhd'];
info = mha_read_header(infpath);
for phase = 1:9
    LMname = [LMpath,'case',num2str(phase),'0.txt'];
    cur_pts = read_POPI_points_file(LMname);
    pts_all_phase(:,:,phase) = round(cur_pts./koef);   
end

min1 = min(min(squeeze(pts_all_phase(:,1,:))));
max1 = max(max(squeeze(pts_all_phase(:,1,:))));
min2 = min(min(squeeze(pts_all_phase(:,2,:))));
max2 = max(max(squeeze(pts_all_phase(:,2,:))));
min3 = min(min(squeeze(pts_all_phase(:,3,:))));
max3 = max(max(squeeze(pts_all_phase(:,3,:))));
%     min_max1 = [ min(pts_mov_max, [], 1)', max(pts_mov_max, [], 1)'];
%     min_max2 = [ min(pts_fix, [], 1)', max(pts_fix, [], 1)'];
%     min_max = [ min(min_max1(:, 1), min_max2(:, 1)), max(min_max1(:, 2), min_max2(:, 2))];
min_max = [min1,max1;min2,max2;min3,max3];

for phase = 0:9
    filepath = 'H:/POPI_model/example/4DCT_MetaImage/';
    init_size = info.Dimensions;
    filename = [filepath,num2str(phase),'0_P.raw'];
    volmov = readrawPOPImeta(filename,init_size);
    volfix=readrawPOPImeta('H:/POPI_model/example/4DCT_MetaImage/10_P.raw',init_size);
    LMpath = 'H:/POPI_model/example/4DLandmarks/';
    LMname = [LMpath,'case',num2str(phase),'0.txt'];
    pts_fix_or = read_POPI_points_file('H:/POPI_model/example/4DLandmarks/case10.txt');
    pts_mov_or = read_POPI_points_file(LMname);
    TIME_e = zeros(2,2,2, 10); 
    filenamesave = [filepath,'case',num2str(phase),'_0.mat'];


spc=[0.976562 0.976562 2];
idx=1;
koef=repmat(spc,[41,1]);
pts_fix=round(pts_fix_or./koef);
pts_mov=round(pts_mov_or./koef);
volmov=volmov+1024;
volfix=volfix+1024;


    volmov = img_thr(volmov, min(volmov(:)), max(volmov(:)), 1);
    volfix = img_thr(volfix, min(volfix(:)), max(volfix(:)), 1);
    
    d = [10, 10, 5];
    crop_v = [ max(1, min_max(1,1) - d(1)), min(size(volmov, 1), min_max(1,2) + d(1)); ... 
               max(1, min_max(2,1) - d(2)), min(size(volmov, 2), min_max(2,2) + d(2)); ... 
               max(1, min_max(3,1) - d(3)), min(size(volmov, 3), min_max(3,2) + d(3));];
%            crop_v = [1,init_size(1);1, init_size(2);1,init_size(3)];
    volmov = crop_data(volmov, crop_v);
    volfix = crop_data(volfix, crop_v);
    

    spc_orig = spc;
    if resize
        bszv = size(volmov);
        spc_tmp = [1, 1, 1];
        volfix = volresize(volfix, round(bszv .* spc .* spc_tmp), 1);
        volmov = volresize(volmov, round(bszv .* spc .* spc_tmp), 1);
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
%     opts.metric = 'loc_cc_fftn_single';
%     opts.metric = 'sad';
    opts.metric_param = [1,1,1] * 2.1;

    opts.scale_metric_param = true;
    opts.isoTV = 0.11;
    opts.csqrt = 5e-3;
    opts.spline_order = 1;
    opts.border_mask = 5;
    opts.max_iters =  80;
    opts.check_gradients = 100*0;
    opts.pix_resolution = spc;

    timer = tic;
    [voldef, Tptv, Kptv] = ptv_register(volmov, volfix, opts);
    TIME_e(use_refinement+1, resize+1, fast_lcc+1, idx) = toc(timer);

    if resize
        Tptv_rsz = cat(4, volresize(Tptv(:,:,:,1), bszv), volresize(Tptv(:,:,:,2), bszv), volresize(Tptv(:,:,:,3), bszv));   
        voldef = volresize(voldef, bszv);
    else
        Tptv_rsz = Tptv;
        
    end
    Tptv = single(Tptv);
    [~, Tptv_rsz] = uncrop_data(voldef, Tptv_rsz, crop_v, init_size);
    save(filenamesave,'Tptv');

    [nr,nc,ns,N] = size(Tptv_rsz);
    xd_kf = single(repmat((1:nr)',1,nc,ns))*spc_orig(1);
    yd_kf = single(repmat((1:nc),nr,1,ns))*spc_orig(2);
    cur = zeros(nr,nc,ns,'single');
    for ly = 1:ns
        cur(:,:,ly) = ones(nr,nc)*ly;
    end
    zd_kf = single(cur)*spc_orig(3);
    griddedInterpolantX = griddedInterpolant(xd_kf,yd_kf,zd_kf,Tptv_rsz(:,:,:,1));
    griddedInterpolantY = griddedInterpolant(xd_kf,yd_kf,zd_kf,Tptv_rsz(:,:,:,2));
    griddedInterpolantZ = griddedInterpolant(xd_kf,yd_kf,zd_kf,Tptv_rsz(:,:,:,3));
    displ(:,1) = griddedInterpolantX(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
    displ(:,2) = griddedInterpolantY(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
    displ(:,3) = griddedInterpolantZ(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
    warpedLMs = displ+pts_fix_or;
    koef = repmat(spc, [size(pts_mov_or, 1), 1]);
    tre = mean(sqrt(sum(((warpedLMs-pts_mov_or)).^2,2)));
    TRE_me(phase) = tre;
    [pt_errs_phys, pts_moved_pix, TRE_phys, TREstd_phys] = DIR_movepoints(pts_mov, pts_fix, Tptv_rsz, spc_orig, []);
    TRE_iso(phase) = mean(TRE_phys);


end

