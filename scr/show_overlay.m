function show_overlay(img_CT, img_GT, mask, index_slice,colormap_name)
% Overlay GT pseudo color image into CT image
% normalization(0-1)
img_CT = img_thr(img_CT, -1024, 0, 1);
max_value_GT = max(img_GT(:));
min_value_GT = min(img_GT(:));
img_GT = img_thr(img_GT, min_value_GT, max_value_GT, 1);

GT_size = size(img_GT);
img_CT_rsz = volresize(img_CT, GT_size, 1);
img_GT_s = img_GT(:,:,index_slice);

img_CT_s = img_CT_rsz(:,:,index_slice);
mask_s = double(mask(:,:,index_slice));
CT_RGB = repmat(img_CT_s,[1,1,3]);
CT_alpha = 1 - mask_s;

img_GT_s_mask = double(img_GT_s).*double(mask_s);
max_value_slice = max(img_GT_s_mask(:));
min_value_slice = min(img_GT_s_mask(:));
img_GT_s_mask = img_thr(img_GT_s_mask, min_value_slice, max_value_slice, 1);

if strcmp(colormap_name,'jet')
    imshow(img_GT_s_mask,'Colormap',jet(255));
else
    if strcmp(colormap_name,'parula')
        imshow(img_GT_s_mask,'Colormap',parula(255));
    else
        if strcmp(colormap_name,'hot')
            imshow(img_GT_s_mask,'Colormap',hot(255));
        else
            if strcmp(colormap_name,'autumn')
                 imshow(img_GT_s_mask,'Colormap',autumn(583));
            end
        end
    end
end
hold on;
h = imshow(CT_RGB,[]);
set(h, 'alphadata', CT_alpha);
