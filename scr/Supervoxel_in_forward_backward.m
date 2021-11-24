function [V_after,C_p_in] = Supervoxel_in_forward_backward(SLIC_Labels_3D,Tptv_rsz,spc,labels_index)
Numlabels = size(labels_index,1);
Tmin_pix = conv_3d_T_from_phys_to_pix(Tptv_rsz, spc);
Tmin_pix_inverse = inverse_field(Tmin_pix);
Tx = Tmin_pix(:,:,:,1);
Ty = Tmin_pix(:,:,:,2);
Tz = Tmin_pix(:,:,:,3);
Tx_inverse = Tmin_pix_inverse(:,:,:,1);
Ty_inverse = Tmin_pix_inverse(:,:,:,2);
Tz_inverse = Tmin_pix_inverse(:,:,:,3);

L_size = size(SLIC_Labels_3D);
index_all = [1:numel(SLIC_Labels_3D)];
[x_index,y_index,z_index] = ind2sub(L_size,index_all);
x_index = x_index';
y_index = y_index';
z_index = z_index';
x_back = round(x_index+Tx_inverse(:));
y_back = round(y_index+Ty_inverse(:));
z_back = round(z_index+Tz_inverse(:));
x_back(find(x_back<1)) = 1;
y_back(find(y_back<1)) = 1;
z_back(find(z_back<1)) = 1;
x_back(find(x_back>L_size(1))) = L_size(1);
y_back(find(y_back>L_size(2))) = L_size(2);
z_back(find(z_back>L_size(3))) = L_size(3);
back_index_all = sub2ind(L_size,x_back,y_back,z_back);

x_forward = round(x_index+Tx(:));
y_forward = round(y_index+Ty(:));
z_forward = round(z_index+Tz(:));
x_forward(find(x_forward <1)) = 1;
y_forward(find(y_forward<1)) = 1;
z_forward(find(z_forward<1)) = 1;
x_forward(find(x_forward >L_size(1))) = L_size(1);
y_forward(find(y_forward>L_size(2))) = L_size(2);
z_forward(find(z_forward>L_size(3))) = L_size(3);
forward_index_all = sub2ind(L_size,x_forward,y_forward,z_forward);

V_after = zeros(Numlabels,1);
C_p_in = cell(Numlabels,1);
for i = 1:Numlabels
    index_cur = find(SLIC_Labels_3D==labels_index(i));
    [Lia,Locb] = ismember(index_cur,back_index_all);
    Locb(find(Locb==0)) = [];
    back_index_cur = Locb;
%     back_index_cur = find(back_index_all==index_cur);
    forward_index_cur = forward_index_all(index_cur);
    sup_in_index = intersect(back_index_cur,forward_index_cur);
    [new_x,new_y,new_z] = ind2sub(L_size,sup_in_index);
    [c_x,c_y,c_z] = ind2sub(L_size,index_cur);
%     Tx_cur = Tx(index_cur);
%     Ty_cur = Ty(index_cur);
%     Tz_cur = Tz(index_cur);
%     new_x = c_x + Tx_cur;
%     new_y = c_y + Ty_cur;
%     new_z = c_z + Tz_cur;
    
    
    
    shp = alphaShape(new_x,new_y,new_z);
    V_after(i) = volume(shp);
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
    area_index = sub2ind(L_size,X_index_cur,Y_index_cur,Z_index_cur);
    if size(unique(c_x),1)==1||size(unique(c_y),1)==1||size(unique(c_z),1)==1
        C_p_in{i} = [];
    else
        C_p_in_cur = double(inShape(shp,X_grid,Y_grid,Z_grid));
        C_p_in_cur = reshape(C_p_in_cur,numel(X_grid),1);
        if isempty(find(C_p_in_cur== 1))
            C_p_in{i} = [];
        else
        
        C_p_in{i} = cat(2,area_index,C_p_in_cur);
        end
    end
end