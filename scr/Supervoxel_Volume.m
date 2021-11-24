% calculate volume
function [V, C_p_ex,supvoxel_sz] = Supervoxel_Volume(SLIC_Labels_3D,labels_index)
Numlabels = size(labels_index,1);
V = zeros(Numlabels,1);
L_size = size(SLIC_Labels_3D);
C_p_ex = cell(Numlabels,1);
supvoxel_sz = zeros(Numlabels,6);
parfor index = 1:Numlabels
    im_index = find(SLIC_Labels_3D == labels_index(index));
    [c_x,c_y,c_z] = ind2sub(L_size,im_index);
    % C_p = cat(2,c_x,c_y,c_z);
    shp = alphaShape(c_x,c_y,c_z);
    V(index) = volume(shp);
    % figure;
    % plot(shp);
    
    [c_x,c_y,c_z] = ind2sub(L_size,im_index);
    
    new_x = c_x ;
    new_y = c_y ;
    new_z = c_z ;
    shp = alphaShape(new_x,new_y,new_z);
    
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
    
    X_index_cur = reshape(X_grid,numel(X_grid),1);
    Y_index_cur = reshape(Y_grid,numel(Y_grid),1);
    Z_index_cur = reshape(Z_grid,numel(Z_grid),1);
    area_index = sub2ind(L_size,X_index_cur,Y_index_cur,Z_index_cur);
    if size(unique(c_x),1)==1||size(unique(c_y),1)==1||size(unique(c_z),1)==1
        C_p_ex{index} = [];
        supvoxel_sz(index,:) = zeros(1,6);
    else
        C_p_ex_cur = double(inShape(shp,X_grid,Y_grid,Z_grid));
        C_p_ex_cur = reshape(C_p_ex_cur,numel(X_grid),1);
        C_p_ex{index} = cat(2,area_index,C_p_ex_cur);
        supvoxel_sz(index,:) =[min_x,max_x,min_y,max_y,min_z,max_z];
    end
    
end
end
