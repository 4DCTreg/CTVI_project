function show_supervoxel_ventilation(vol_ex,vent_img_sup_cur)
Slice = round(size(vent_img_sup_cur,2)/2);
img_2D = squeeze(vol_ex(:,Slice,:));
img_2D = imrotate(img_2D,90);
img_2D = img_thr(img_2D,min(img_2D(:)),max(img_2D(:)),1);
vent_img_2D = squeeze(vent_img_sup_cur(:,Slice,:));
vent_img_2D = imrotate(vent_img_2D,90);

mask = vent_img_2D;
mask(find(mask~=0))=1;
vent_index = reshape(vent_img_2D,[numel(vent_img_2D),1]);
vent_index_sort = sort(vent_index);
min_p_value = vent_index_sort(round(0.005*numel(vent_img_2D)));
max_p_value = vent_index_sort(round(0.995*numel(vent_img_2D)));
vent_img_2D = img_thr(vent_img_2D,min_p_value,max_p_value,1);
vent_img_2D = vent_img_2D.*mask;
CT_alpha = 1-mask;
k1 = unique(vent_img_2D);
dx = 1/length(k1);
inter = [0:dx:1];

for un_idx = 1:length(k1)
    vent_img_2D(find(vent_img_2D==k1(un_idx))) = inter(un_idx);
end
mycolormap = hot(length(k1));

figure


imagesc(vent_img_2D);
colormap(mycolormap);
hold on;
img_2D_RGB = repmat(img_2D,[1,1,3]);
h = imshow(img_2D_RGB,[]);
set(h,'alphadata',CT_alpha);