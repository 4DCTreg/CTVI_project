clear;
addpath(genpath('IOfiles/'));
volfix_mask=readrawPOPI('H:\POPI_model\example\4DMask_metalmage\60-air-body-lungs.raw');
volfix_mask(volfix_mask<2)=0;
volfix_mask(volfix_mask==2)=1;
volinspiration=readrawPOPImeta('H:\POPI_model\example\4DCT_MetaImage\10_P.raw');
volexpiration=readrawPOPImeta('H:\POPI_model\example\4DCT_MetaImage\60_P.raw');
volfix=volinspiration;
volmov=volexpiration;
conoral_or=volfix(:,:,70);
alpha_mask=volfix_mask;
alpha_mask(alpha_mask==0)=2;
alpha_mask(alpha_mask==1)=0;
alpha_mask(alpha_mask==2)=1;
coronal_mask=alpha_mask(:,:,70);

max_val=max(conoral_or(:));
min_val=min(conoral_or(:));
coronal_color=zeros(482,360,3);
coronal_color(:,:,1)=conoral_or;
coronal_color(:,:,2)=conoral_or;
coronal_color(:,:,3)=conoral_or;
coronal_color=(coronal_color-min_val)/(max_val-min_val);
% imshow(coronal_color);
% figure
% imshow(conoral_or,[]);
% % set(h,'AlphaData',coronal_mask);
% hold on
meanfilter=(1/27)*ones(3,3,3);
volfix_mean=imfilter(volfix,meanfilter,'symmetric');
% volfix_mean=padarray(volfix,[1,1,1],'symmetric');
load('dvf.mat');
[r,c,s]=size(volmov);
spc=[0.976562 0.976562 2];
[Y,X,Z]=meshgrid(0:1:c-1,0:1:r-1,0:1:s-1);
X=X.*spc(1)+spc(1)/2;
Y=Y.*spc(2)+spc(2)/2;
Z=Z.*spc(3)+spc(3)/2;

mov_x=X+Tptv_rsz(:,:,:,1);
mov_y=Y+Tptv_rsz(:,:,:,2);
mov_z=Z+Tptv_rsz(:,:,:,3);

x_index=round((mov_x-spc(1)/2)/spc(1))+1;
y_index=round((mov_y-spc(2)/2)/spc(2))+1;
z_index=round((mov_z-spc(3)/2)/spc(3))+1;

mov_num_index_matric=ones(r,c,s);
mov_num_index=zeros(r,c,s);
x_index(x_index<=0)=1;
x_index(x_index>r)=r;
y_index(y_index<=0)=1;
y_index(y_index>c)=c;
z_index(z_index<=0)=1;
z_index(z_index>s)=s;
inhale_average=zeros(r,c,s);
% node_index=zeros(r*c*s,6);
% for index=1:r*c*s
%     [i,j,k]=ind2sub([r,c,s],index);
%     node_index(index,1)=i;
%     node_index(index,2)=j;
%     node_index(index,3)=k;
%     node_index(index,4)=x_index(i,j,k);
%     node_index(index,5)=y_index(i,j,k);
%     node_index(index,6)=z_index(i,j,k);
%             
% end
% node_ex=node_index(:,4:6);
% node_ins=node_index(:,1:3);
for i=1:r
    for j=1:c
        for k=1:s
            x=x_index(i,j,k);
            y=y_index(i,j,k);
            z=z_index(i,j,k);
            inhale_average(x,y,z)=volinspiration(i,j,k)+inhale_average(x,y,z);
            mov_num_index(x,y,z)=mov_num_index(x,y,z)+1;
%             if volfix_mask(i,j,k)==1
%             cur_ex=[i,j,k];
%             index=ismember(node_ex,cur_ex,'rows');
%             ins_all=find(index==1);
%             
%             if ins_all>0
%                 n=size(ins_all,1);
%                 for num=1:n
%                     index_ins=ins_all(num);
%                     inhale_average(i,j,k)=volinspiration(node_ins(index_ins,1),node_ins(index_ins,2),node_ins(index_ins,3))+inhale_average(i,j,k);
%                 end
%                 inhale_average(i,j,k)=inhale_average(i,j,k)/n;
%                 mov_num_index(i,j,k)=n;
%             else
%                 inhale_average(i,j,k)=volexpiration(i,j,k);
%                 mov_num_index(i,j,k)=0;
%             end
            
%             x=x_index(i,j,k);
%             y=y_index(i,j,k);
%             z=z_index(i,j,k);
%             inhale_average(i,j,k)=volmov(x,y,z)+inhale_average(i,j,k);
%             
%             mov_num_index(i,j,k)=mov_num_index_matric(x,y,z)+mov_num_index(i,j,k);
            
%             end
                
        end
    end
end

index_non=find(mov_num_index==0);
inhale_average(index_non)=volexpiration(index_non);
mov_num_index(index_non)=1;
inhale_mean=inhale_average./mov_num_index;

% inhale_mean=imfilter(inhale_average,meanfilter,'symmetric');
inhale_mean_mask=inhale_mean.*volfix_mask;
inhale_mean_mask(inhale_mean_mask==-1000)=-1001;
mask=volfix_mask;
mask(mask==0)=0.5;
volfix_mean_mask=volfix_mean.*mask;
volfix_mean_mask(volfix_mean_mask==0)=1;
% inhale_average(inhale_average==-1000)=-1001;
% volfix(volfix==0)=0.1;
volfix_mean(volfix_mean==0)=0.01;

ventilation_image_nomask=1000*(inhale_mean_mask-volfix_mean_mask)./(volfix_mean_mask.*(1000+inhale_mean_mask));
% ventilation_image_nomask=1000*(inhale_mean-volfix_mean)./(volfix_mean.*(1000+inhale_mean));
ventilation_image=imfilter(ventilation_image_nomask,meanfilter,'symmetric');
ventilation_image_mask=ventilation_image.*volfix_mask;
ventilation_image_mask(ventilation_image_mask>5)=5;
ventilation_image_mask(ventilation_image_mask<-5)=-5;
coronal=ventilation_image_mask(:,:,70);
coronal=(coronal+5);
coronal_mask_or=volfix_mask(:,:,70);
map=colormap(jet(256));
space=10/255;
coronal_color_lung=zeros(482,360,3);
exhale_mask=volfix_mask(:,:,70);
for i=1:482
    for j=1:360
        if exhale_mask(i,j)==1           
            index=fix(coronal(i,j)/space)+1;
            coronal_color_lung(i,j,1)=map(index,1);
            coronal_color_lung(i,j,2)=map(index,2);
            coronal_color_lung(i,j,3)=map(index,3);    
        else            
            coronal_color_lung(i,j,1)=coronal_color(i,j,1);
            coronal_color_lung(i,j,2)=coronal_color(i,j,2);
            coronal_color_lung(i,j,3)=coronal_color(i,j,3);
        end
    end
end
imshow(coronal_color_lung);
% x=imshow(coronal,map);
% h=imagesc(coronal);
% colormap jet
% set(h,'AlphaData',coronal_mask_or);
colorbar;
% colormap jet
% ventilation_image=1000*(inhale_mean_mask-volfix_mean_mask)./(volfix_mean_mask.*(1000+inhale_mean_mask));
% 
% ventilation_image(isnan(ventilation_image))=0;
% ventilation_image(ventilation_image==inf)=0;
% ventilation_image(ventilation_image==-inf)=0;
