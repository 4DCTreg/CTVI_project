clear;
load('HU_debug.mat');
mask_index = find(vol_ex<=-250);
no_mask = find(vol_ex>-250);
Tmin_pix = conv_3d_T_from_phys_to_pix(T_rsz, spc);
Tmin_pix_inverse = inverse_field(Tmin_pix);
vol_size = size(vol_ex);
Tx = T_rsz(:,:,:,1);
Ty = T_rsz(:,:,:,2);
Tz = T_rsz(:,:,:,3);
Tx_inverse = Tmin_pix_inverse(:,:,:,1);
Ty_inverse = Tmin_pix_inverse(:,:,:,2);
Tz_inverse = Tmin_pix_inverse(:,:,:,3);
count = zeros(vol_size);
avg_img_in = zeros(vol_size);
for idx = 1:vol_size(1)*vol_size(2)*vol_size(3)
    [i,j,k] = ind2sub(vol_size,idx);
    dx = Tx_inverse(idx);
    dy = Ty_inverse(idx);
    dz = Tz_inverse(idx);
    new_i = round(i+dx);
    if new_i<=0
        new_i = 1;
    end
    if new_i>vol_size(1)
        new_i = vol_size(1);
    end
    new_j = round(j+dy);
    if new_j<=0
        new_j = 1;
    end
    if new_j>vol_size(2)
        new_j = vol_size(2);
    end
    new_k = round(k+dz);
    if new_k<=0
        new_k = 1;
    end
    if new_k>vol_size(3)
        new_k = vol_size(3);
    end
    new_idx = sub2ind(vol_size,new_i,new_j,new_k);
    avg_img_in(new_idx) = avg_img_in(new_idx)+vol_in(idx);
    count(new_idx) = count(new_idx)+1;
end
no_map = find(count==0);
[no_x,no_y,no_z] = ind2sub(vol_size,no_map);
no_map_dx = Tx(no_map);
no_map_dy = Ty(no_map);
no_map_dz = Tz(no_map);
map_i = no_x+round(no_map_dx);
map_j = no_y+round(no_map_dy);
map_k = no_z+round(no_map_dz);
map_i(find(map_i<=0)) = 1;
map_j(find(map_j<=0)) = 1;
map_k(find(map_k<=0)) = 1;
map_i(find(map_i>vol_size(1))) = vol_size(1);
map_j(find(map_j>vol_size(2))) = vol_size(2);
map_k(find(map_k>vol_size(3))) = vol_size(3);
in_idx = sub2ind(vol_size,map_i,map_j,map_k);
avg_img_in(no_map) = vol_in(in_idx);
count(find(count==0)) = 1;
vol_avg = avg_img_in./count;

vent_img = 1000*(vol_avg-vol_ex)./(vol_ex.*(1000+vol_avg));
vent_img_mask = vent_img;
vent_img_mask(no_mask) = 0;
vol_img_mask = img_thr(vol_img_mask, min(vol_img_mask(:)), max(vol_img_mask(:)), 100);