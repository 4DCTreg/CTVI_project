clear
addpath(genpath('../mutils/My/'));
addpath(genpath('../ptv'));
basepth = 'E:/data_prj/dir_dataset/DIR_files/';
addpath(genpath('G:/SIFT3D/SIFT3D 1.4.5/'));
addpath(genpath('../../../scr/'));


TIME_e = zeros(2,2,2, 10); 
TRE_e = zeros(2,2,2, 10);
TRE_SIFT = zeros(10,1);
TRE_SIFT_or = zeros(10,1);

use_refinement = 0 ;
resize = 1;
fast_lcc =  1;


for idx = 1:10
    pts_struct = DIR_get_all_points_for_the_case(idx, basepth);
    [volmov, spc] = read_DIR_volume_4dCT(idx, 5, basepth);
    volmov = double(volmov);
    pts_mov = pts_struct.extreme.e;
    [volfix, spc] = read_DIR_volume_4dCT(idx, 0, basepth);
    volfix = double(volfix);
    pts_fix = pts_struct.extreme.b;

    volmov = img_thr(volmov, 80, 900, 1);
    volfix = img_thr(volfix, 80, 900, 1);
    
    % crop images 
    init_size = size(volmov);
    min_max1 = [ min(pts_mov, [], 1)', max(pts_mov, [], 1)'];
    min_max2 = [ min(pts_fix, [], 1)', max(pts_fix, [], 1)'];
    min_max = [ min(min_max1(:, 1), min_max2(:, 1)), max(min_max1(:, 2), min_max2(:, 2))];
    d = [10, 10, 5];
    crop_v = [ max(1, min_max(1,1) - d(1)), min(size(volmov, 1), min_max(1,2) + d(1)); ... 
               max(1, min_max(2,1) - d(2)), min(size(volmov, 2), min_max(2,2) + d(2)); ... 
               max(1, min_max(3,1) - d(3)), min(size(volmov, 3), min_max(3,2) + d(3));];
    volmov = crop_data(volmov, crop_v);
    volfix = crop_data(volfix, crop_v);
    units = spc;
    % Detect keypoints
    keys_in = detectSift3D(volmov, 'units', units);
    keys_ex = detectSift3D(volfix, 'units', units);
    % Extract descriptors
    [desc_in, coords_in] = extractSift3D(keys_in);
    [desc_ex, coords_ex] = extractSift3D(keys_ex);
    LM = coords_in;
    LF = coords_ex;
    correspondences = vl_ubcmatch(desc_in',desc_ex');
    [keyM,keyF] = calc_affine(LM,LF,correspondences);
    pts_mov_SIFT = keyM;
    pts_mov_SIFT(:,1) = pts_mov_SIFT(:,1) + crop_v(1,1) - 1;
    pts_mov_SIFT(:,2) = pts_mov_SIFT(:,2) + crop_v(2,1) - 1;
    pts_mov_SIFT(:,3) = pts_mov_SIFT(:,3) + crop_v(3,1) - 1;
    pts_fix_SIFT = keyF;
    pts_fix_SIFT(:,1) = pts_fix_SIFT(:,1) + crop_v(1,1) - 1;
    pts_fix_SIFT(:,2) = pts_fix_SIFT(:,2) + crop_v(2,1) - 1;
    pts_fix_SIFT(:,3) = pts_fix_SIFT(:,3) + crop_v(3,1) - 1;
    
    koef = repmat(spc, [size(pts_mov_SIFT, 1), 1]);
    pt_errs_phys = sqrt( sum((  (pts_mov_SIFT - pts_fix_SIFT).*koef  ).^2, 2) );
%     index_1 = find(pt_errs_phys==0);
%     index_2 = find(pt_errs_phys > 15);
%     index = [index_1;index_2];
%     pts_fix_SIFT(index,:) = [];
%     pts_mov_SIFT(index,:) = [];
%     pt_errs_phys(index) = [];
    index_all = exclude_neighbor(pts_fix_SIFT,pts_mov_SIFT,spc);
    pts_fix_SIFT(index_all,:) = [];
    pts_mov_SIFT(index_all,:) = [];
    pt_errs_phys(index_all) = [];
    pt_errs_phys_SIFT_or = pt_errs_phys;
    TRE_phys_SIFT = mean(pt_errs_phys);
    TRE_SIFT_or(idx) = TRE_phys_SIFT;
%     TREstd_phys_SIFT = std(sqrt( sum((  (pts_mov_SIFT - pts_fix_SIFT).*koef  ).^2, 2) ));

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
    [~, Tptv_rsz] = uncrop_data(voldef, Tptv_rsz, crop_v, init_size);
    % move points and measure TRE
    [pt_errs_phys, pts_moved_pix, TRE_phys, TREstd_phys] = DIR_movepoints(pts_mov, pts_fix, Tptv_rsz, spc_orig, []);
    TREs(idx) = mean(TRE_phys);
    
    [pt_errs_phys_or, TRE_phys_or, TREstd_phys_or] = DIR_or(pts_mov, pts_fix, spc_orig, []);
    
    
    [pt_errs_phys_SIFT, pts_moved_pix, TRE_phys, TREstd_phys] = DIR_movepoints(pts_mov_SIFT, pts_fix_SIFT, Tptv_rsz, spc_orig, []);
%     pt_errs_phys_SIFT = pt_errs_phys;
    TRE_SIFT(idx) = mean(TRE_phys);

    TRE_e(use_refinement+1, resize+1, fast_lcc+1, idx) = mean(TRE_phys);

end
