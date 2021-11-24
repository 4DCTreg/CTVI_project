% This code is used to visualization of 3D deformation field in kalman
% estimate flow
clear
filePathiso_rsz = 'H:/continus_reg/ObservedData/isoPTV_continus_rsz/';
method = 'isoPTV_norsz';
idx = 5;
T_ob = [];   
for phase = 5
    fieldFilename = [filePathiso_rsz,'case',num2str(idx),'/','T',num2str(phase),'_0','.mat'];
    load(fieldFilename);
    T_ob(:,:,:,:) = Tptv;
end
init_sz = size(T_ob);
m = 7; n = 7; h = 7; 
bszv = [m,n,h];
T_iso = cat(4, volresize(T_ob(:,:,:,1), bszv), ...
       volresize(T_ob(:,:,:,2), bszv), volresize(T_ob(:,:,:,3), bszv)); 

[X_grid,Y_grid,Z_grid] = ndgrid(1:m,1:n,1:h) ;
% map = [1 0.72549 0.05882];
% map = [0.63529	0.80392	0.35294];
% map = [0.3098	0.58039	0.80392];
% map = [0.93333 0 0]; 
map = [0.09412 0.4549 0.80392];
map = [0 0 0];

colormap(map);
cm = ones(7,7);
for lz = 1:7
    Z_grid(:,:,lz) = ones(7,7)*(7-lz+1);
end
Dx = X_grid+T_iso(:,:,:,1)./(init_sz(1)/7);
Dy = Y_grid+T_iso(:,:,:,2)./(init_sz(2)/7);
Dz = Z_grid+T_iso(:,:,:,3)./(init_sz(3)/7);
for lz = 1:7
    mesh(Dx(:,:,lz),Dy(:,:,lz),Dz(:,:,lz),cm,'LineWidth',8);
    hold on
end

for lx = 1:7

    mesh(squeeze(Dx(lx,:,:)),squeeze(Dy(lx,:,:)),squeeze(Dz(lx,:,:)),cm,'LineWidth',2.5);
    hold on
end
for ly = 1:7

    mesh(squeeze(Dx(:,ly,:)),squeeze(Dy(:,ly,:)),squeeze(Dz(:,ly,:)),cm,'LineWidth',2.5);
    hold on
end
alpha(0.3);
axis off
view(-56,27);
