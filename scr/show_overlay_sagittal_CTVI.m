function show_overlay_sagittal_CTVI(img_CT, img_GT, info_GT, mask, index_slice,colormap_name)
% Overlay GT pseudo color image into CT image
% normalization(0-1)
% sz_gt = info_GT.Dimensions;
img_CT = double(img_CT);
img_GT = double(img_GT);
mask = double(mask);
units = info_GT.PixelDimensions;

bszv = info_GT.Dimensions;
spc_tmp = [1, 1, 1];
% img_GT = volresize(img_GT, round(bszv .* units .* spc_tmp), 1);
img_GT = imresize3(img_GT, round(bszv .* units .* spc_tmp), 'linear');
mask = imresize3(mask, round(bszv .* units .* spc_tmp), 'nearest');
%         vol_ex_nonor = volresize(vol_ex_nonor, round(bszv .* units .* spc_tmp), 1);
%         vol_in_nonor = volresi(vol_in_nonor, round(bszv .* units .* spc_tmp), 1);
%         avg_img = volresize(avg_img, round(bszv .* units .* spc_tmp), 1);
%         spc = [1,1,1] ./ spc_tmp;
new_sz = size(img_GT);
img_GT_index = reshape(img_GT,[new_sz(1)*new_sz(2)*new_sz(3),1]);
img_GT_sort = sort(img_GT_index);
p_value = round(img_GT_sort(round(0.99*new_sz(1)*new_sz(2)*new_sz(3))));
p_value_min = round(img_GT_sort(round(0.02*new_sz(1)*new_sz(2)*new_sz(3))));
img_CT = img_thr(img_CT, -1024, 0, 1);
max_value_GT = max(img_GT(:));
min_value_GT = min(img_GT(:));
img_GT = img_thr(img_GT, min_value_GT, max_value_GT, 1);

GT_size = size(img_GT);
img_CT_rsz = volresize(img_CT, GT_size, 1);
img_GT_s = squeeze(img_GT(:,index_slice,:));
img_GT_s = imrotate(img_GT_s,90);




img_CT_s = squeeze(img_CT_rsz(:,index_slice,:));
img_CT_s = imrotate(img_CT_s,90);
mask_s = squeeze(double(mask(:,index_slice,:)));
mask_s = imrotate(mask_s,90);
CT_RGB = repmat(img_CT_s,[1,1,3]);
CT_alpha = 1 - mask_s;

img_GT_s_mask = double(img_GT_s).*double(mask_s);
max_value_slice = max(img_GT_s_mask(:));
min_value_slice = min(img_GT_s_mask(:));
img_GT_s_mask = img_thr(img_GT_s_mask, min_value_slice, max_value_slice, 1);
u_value = unique(img_GT_s_mask);
dx = 1/length(u_value);
inter = [0:dx:1];
for un_idx = 1:length(u_value)
    img_GT_s_mask(find(img_GT_s_mask==u_value(un_idx))) = inter(un_idx);
end

if strcmp(colormap_name,'jet')
    imshow(img_GT_s_mask,'Colormap',jet(255));
else
    if strcmp(colormap_name,'parula')
        imshow(img_GT_s_mask,'Colormap',parula(255));
    else
        if strcmp(colormap_name,'hot')
            imshow(img_GT_s_mask,'Colormap',hot(1000));
        else
            if strcmp(colormap_name,'autumn')
                imshow(img_GT_s_mask,'Colormap',autumn(size(u_value,1)));
            end
        end
    end
end
hold on;
h = imshow(CT_RGB,[]);
set(h, 'alphadata', CT_alpha);
end
