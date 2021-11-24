function [pts_in,pts_ex] = SIFT_extract_landmarks(vol_in, vol_ex, units, crop_v,avg_mask)
spc = units;
    keys_in = detectSift3D(vol_in, 'units', units);
    keys_ex = detectSift3D(vol_ex, 'units', units);
    % Extract descriptors
    [desc_in, coords_in] = extractSift3D(keys_in);
    [desc_ex, coords_ex] = extractSift3D(keys_ex);
    L_in = coords_in;
    L_ex = coords_ex;
    correspondences = vl_ubcmatch(desc_in',desc_ex',2);
    [key_in,key_ex] = calc_affine(L_in,L_ex,correspondences);
    pts_in = key_in;
    pts_in(:,1) = pts_in(:,1) + crop_v(1,1) - 1;
    pts_in(:,2) = pts_in(:,2) + crop_v(2,1) - 1;
    pts_in(:,3) = pts_in(:,3) + crop_v(3,1) - 1;
    pts_ex = key_ex;
    pts_ex(:,1) = pts_ex(:,1) + crop_v(1,1) - 1;
    pts_ex(:,2) = pts_ex(:,2) + crop_v(2,1) - 1;
    pts_ex(:,3) = pts_ex(:,3) + crop_v(3,1) - 1;
    
    koef = repmat(spc, [size(pts_in, 1), 1]);
    pt_errs_phys = sqrt( sum((  (pts_in - pts_ex).*koef  ).^2, 2) );
    index_1 = find(pt_errs_phys==0);
    index_2 = find(pt_errs_phys > 15);
    index = [index_1;index_2];
    index = index_1;
    pts_ex(index,:) = [];
    pts_in(index,:) = [];
    pt_errs_phys(index) = [];

    index_all = exclude_neighbor(pts_ex,pts_in,spc);
    pts_ex(index_all,:) = [];
    pts_in(index_all,:) = [];
    pt_errs_phys(index_all) = [];
    index_mask = [];
    for idx = 1:size(pts_ex,1)
        mask_cur = avg_mask(pts_ex(idx,1),pts_ex(idx,2),pts_ex(idx,3));
        if mask_cur == 0
            index_mask = [index_mask;idx];
        end
    end
    pts_ex(index_mask,:) = [];
    pts_in(index_mask,:) = [];


end
