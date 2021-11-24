function [V_before,V_after,N_in_after,C_p_in] = Supervoxel_inhale(SLIC_Labels_3D,Tptv_rsz,spc,labels_index)
Numlabels = size(labels_index,1);
Tmin_pix = conv_3d_T_from_phys_to_pix(Tptv_rsz, spc);
Tx = Tmin_pix(:,:,:,1);
Ty = Tmin_pix(:,:,:,2);
Tz = Tmin_pix(:,:,:,3);
L_size = size(SLIC_Labels_3D);
V_after = zeros(Numlabels,1);
V_before = zeros(Numlabels,1);
C_p_in = cell(Numlabels,1);
N_in_after = zeros(Numlabels,1);


Tmin_pix_inverse = inverse_field(Tmin_pix);
vol_size = L_size;
% Tx = T_rsz(:,:,:,1);
% Ty = T_rsz(:,:,:,2);
% Tz = T_rsz(:,:,:,3);
Tx_inverse = Tmin_pix_inverse(:,:,:,1);
Ty_inverse = Tmin_pix_inverse(:,:,:,2);
Tz_inverse = Tmin_pix_inverse(:,:,:,3);


vol_in_to_ex = zeros(vol_size);
parfor idx = 1:vol_size(1)*vol_size(2)*vol_size(3)
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
vol_in_to_ex(idx) =new_idx; 

end




for i = 1:Numlabels
    sup_in_all = [];
    index_cur = find(SLIC_Labels_3D==labels_index(i));
    l_m = ismember(vol_in_to_ex,index_cur);
    sup_in_all = find(l_m==1);
%     for indx = 1:numel(index_cur)
%     sup_inverse = find(vol_in_to_ex==index_cur(indx));
%     sup_in_all = [sup_in_all;sup_inverse];
%     end
    
    V_before(i) = numel(index_cur);
    [c_x,c_y,c_z] = ind2sub(L_size,index_cur);
    Tx_cur = Tx(index_cur);
    Ty_cur = Ty(index_cur);
    Tz_cur = Tz(index_cur);
    new_x = round(c_x + Tx_cur);
    new_y = round(c_y + Ty_cur);
    new_z = round(c_z + Tz_cur);
    
    t_x_min = new_x<=0;
    new_x(t_x_min) = 1;
    t_x_max = new_x>vol_size(1);
    new_x(t_x_max) = vol_size(1);
    
    t_y_min = new_y<=0;
    new_y(t_y_min) = 1;
    t_y_max = new_y>vol_size(2);
    new_y(t_y_max) = vol_size(2);
    
    t_z_min = new_z<=0;
    new_z(t_z_min) = 1;
    t_z_max = new_z>vol_size(3);
    new_z(t_z_max) = vol_size(3);
   
    
    
    area_index = sub2ind(L_size,new_x,new_y,new_z);
%     shp = alphaShape(new_x,new_y,new_z);
%     V_after(i) = volume(shp);
    %     figure;
    %     plot(shp);
    min_x = floor(min(new_x));
    if min_x <= 0
        min_x = 1;
    end
    max_x = ceil(max(new_x));
    if max_x > L_size(1)
        max_x = L_size(1);
    end
    
    
    min_y = floor(min(new_y));
    if min_y <= 0
        min_y = 1;
    end
    
    max_y = ceil(max(new_y));
    if max_y > L_size(2)
        max_y = L_size(2);
    end
    min_z = floor(min(new_z));
    if min_z<=0
        min_z = 1;
    end
    max_z = ceil(max(new_z));
    if max_z > L_size(3)
        max_z = L_size(3);
    end
    [Y_grid,X_grid,Z_grid] = meshgrid(min_y:max_y,min_x:max_x,min_z:max_z);
    sp_size = size(X_grid);
    X_index_cur = reshape(X_grid,numel(X_grid),1);
    Y_index_cur = reshape(Y_grid,numel(Y_grid),1);
    Z_index_cur = reshape(Z_grid,numel(Z_grid),1);
%     area_index = sub2ind(L_size,X_index_cur,Y_index_cur,Z_index_cur);
    if size(unique(c_x),1)==1||size(unique(c_y),1)==1||size(unique(c_z),1)==1
        C_p_in{i} = sup_in_all;
    else
%         C_p_in_cur = double(inShape(shp,X_grid,Y_grid,Z_grid));
%         C_p_in_cur = reshape(C_p_in_cur,numel(X_grid),1);
%         N_in_after(i) = numel(find(C_p_in_cur==1));
        sup_in_all = [sup_in_all;area_index];
        C_p_in{i} = unique(sup_in_all);
    end
    V_after(i)= numel(C_p_in{i});
end