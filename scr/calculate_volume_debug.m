% calculate volume
index = find(SLIC_Labels_3D==900);
L_size = size(SLIC_Labels_3D);
[c_x,c_y,c_z] = ind2sub(L_size,index);
C_p = cat(2,c_x,c_y,c_z);
shp = alphaShape(c_x,c_y,c_z);
V = volume(shp);
plot(shp);
