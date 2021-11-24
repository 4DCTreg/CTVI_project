load('debug_HU.mat');
Numlabels = size(labels_index,1);
im_sz = size(vol_ex);
Tmin_pix = conv_3d_T_from_phys_to_pix(Tptv_rsz, spc);
Tmin_pix_inverse = inverse_field(Tmin_pix);
Tx = Tmin_pix(:,:,:,1);
Ty = Tmin_pix(:,:,:,2);
Tz = Tmin_pix(:,:,:,3);
Tx_inverse = Tmin_pix_inverse(:,:,:,1);
Ty_inverse = Tmin_pix_inverse(:,:,:,2);
Tz_inverse = Tmin_pix_inverse(:,:,:,3);

ex_mask_1 = find(vol_ex>0);
ex_mask_2 = find(vol_ex<=-1000);
vol_ex_sup = vol_ex;
vol_ex_sup(ex_mask_1) = 0;
vol_ex_sup(ex_mask_2) = -1000;
in_mask_1 = find(vol_in>0);
in_mask_2 = find(vol_in<=-1000);
vol_in_sup = vol_in;
vol_in_sup(in_mask_1) = 0;
vol_in_sup(in_mask_2) = -1000;


vent_img = zeros(im_sz);
count = zeros(im_sz);
avg_img_in = zeros(im_sz);
for idx = 1:numel(vol_in)
    [i_in,j_in,k_in] = ind2sub(im_sz,idx);
    dx = Tx_inverse(idx);
    dy = Ty_inverse(idx);
    dz = Tz_inverse(idx);
    new_i = round(i_in+dx);
    if new_i<=0
        new_i = 1;
    end
    if new_i>im_sz(1)
        new_i = im_sz(1);
    end
    new_j = round(j_in+dy);
    if new_j<=0
        new_j = 1;
    end
    if new_j>im_sz(2)
        new_j = im_sz(2);
    end
    new_k = round(k_in+dz);
    if new_k<=0
        new_k = 1;
    end
    if new_k>im_sz(3)
        new_k = im_sz(3);
    end
    new_idx = sub2ind(im_sz,new_i,new_j,new_k);
    avg_img_in(new_idx) = avg_img_in(new_idx)+vol_in_sup(idx);
    count(new_idx) = count(new_idx)+1;
end
no_map = find(count==0);
[no_x,no_y,no_z] = ind2sub(im_sz,no_map);
no_map_dx = Tx(no_map);
no_map_dy = Ty(no_map);
no_map_dz = Tz(no_map);
map_i = no_x+round(no_map_dx);
map_j = no_y+round(no_map_dy);
map_k = no_z+round(no_map_dz);
map_i(find(map_i<=0)) = 1;
map_j(find(map_j<=0)) = 1;
map_k(find(map_k<=0)) = 1;
map_i(find(map_i>im_sz(1))) = im_sz(1);
map_j(find(map_j>im_sz(2))) = im_sz(2);
map_k(find(map_k>im_sz(3))) = im_sz(3);
in_idx = sub2ind(im_sz,map_i,map_j,map_k);
avg_img_in(no_map) = vol_in_sup(in_idx);
count(find(count==0)) = 1;
vol_avg = avg_img_in./count;


for id = 1:Numlabels
    vol_ex_index = find(SLIC_Labels_3D==labels_index(id));
    sup_num = size(vol_ex_index,1);
    vol_ex_HU_sup = vol_ex_sup(vol_ex_index);
    vol_ex_avg_cur = mean(vol_ex_HU_sup);
    vol_in_HU_sup = vol_avg(vol_ex_index);
    vol_in_avg_cur = mean(vol_in_HU_sup);
    denominator = 1000+vol_in_avg_cur;
    if denominator==0
        denominator = 1;
    end
    vent_sup = 1000*(vol_in_avg_cur-vol_ex_avg_cur)/denominator;
    vent_img(vol_ex_index) = vent_sup;
    %     volin_sup = C_p_in{id};
    %     C_p_ex_cur = C_p_ex{id};
    %     if isempty(volin_sup)||isempty(C_p_ex_cur)
    %         continue
    %     else
    %         volin_index = volin_sup(:,1);
    %         volin_index_mask = volin_sup(:,2);
    %         volex_air = (1000-vol_ex_sup(vol_ex_index));
    %         volex_air_sum = sum(volex_air(:));
    %         volin_air = ((1000-vol_in_sup(volin_index))).*volin_index_mask;
    %         volin_air_sum = sum(volin_air(:));
    %         sup_vent = volin_air_sum - volex_air_sum;
    %
    %         n_voxel = size(find(C_p_ex_cur(:,2)==1),1);
    %         vent_voxel = zeros(n_voxel,1);
    %         sup_jac_coef = vent_jac_interp(C_p_ex_cur(:,1));
    %         sup_no_mask = find(C_p_ex_cur(:,2)==0);
    %         sup_jac_coef(sup_no_mask)=[];
    %         vol_sup_ex_index = C_p_ex_cur(:,1);
    %         vol_sup_ex_index(sup_no_mask) = [];
    %         %     sup_coef_min = min(sup_jac_coef);
    %         %     sup_coef_max = max(sup_jac_coef);
    %         %     sup_jac_coef = img_thr(sup_jac_coef, sup_coef_min, sup_coef_max, 5);
    %         sp_sz = supvoxel_sz_ex(id,:);
    %         objf = @(X)tv_grad(X, sup_vent, sup_jac_coef,C_p_ex_cur,sp_sz,tv_reg);
    %         options = [];
    %         options.display = 'off';
    %         % options.LS_init = 4;
    %         options.MaxIter = 80;
    %         options.MaxFunEvals = round(80 * 1.3);
    %         options.Method = 'lbfgs';
    %
    %         options.DerivativeCheck = 'off';
    %         options.Damped = 0;
    %         options.LS_type = 0; options.LS_init = 8; %
    %         %     options.LS_type = 0; options.LS_init = 2; %
    %         %     options.LS_type = 1; options.LS_init = 2; %
    %         options.optTol = 1e-8;
    %         options.progTol = 1e-4;
    %
    %         options.Corr = 85;
    %         if numel(vent_voxel) > (100^3)
    %             options.Corr = 30;
    %             %         options.Corr
    %         end
    %         %     options.progTol = 1e-30;
    %         %     options.optTol = 1e-30;
    %
    %         options.A_step_length = 1/2.4;
    %         [xmin, fmin, ~, outp] = minFuncJoker(objf, vent_voxel(:), options);
    %         if sum(xmin)==0
    %             ncoef = 1;
    %         else
    %             ncoef = sup_vent/sum(xmin);
    %         end
    %         xmin = xmin*ncoef;
    %
    %         vent_img(vol_sup_ex_index) = xmin;
    %     end
end
vent_img_mask = vent_img;
vent_index = reshape(vent_img_mask,[im_sz(1)*im_sz(2)*im_sz(3),1]);
vent_index_sort = sort(vent_index);
min_p_value = vent_index_sort(round(0.01*im_sz(1)*im_sz(2)*im_sz(3)));
max_p_value = vent_index_sort(round(0.95*im_sz(1)*im_sz(2)*im_sz(3)));
vent_img_mask(find(vent_img_mask < min_p_value)) = min_p_value;
vent_img_mask(find(vent_img_mask > max_p_value)) = max_p_value;


vent_img_mask = img_thr(vent_img_mask, min(vent_img_mask(:)), max(vent_img_mask(:)), 1500);
vent_img = vent_img_mask;

