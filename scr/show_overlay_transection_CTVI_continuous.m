function show_overlay_transection_CTVI_continuous(img_CT, img_GT, info_GT, mask, index_slice,colormap_name,num_phase)
% Overlay GT pseudo color image into CT image
% normalization(0-1)
% sz_gt = info_GT.Dimensions;
min_p_value = min(img_GT(:));
max_p_value = max(img_GT(:));
img_GT = img_thr(img_GT, min_p_value, max_p_value, 1500);
units = info_GT.PixelDimensions;

bszv = info_GT.Dimensions;
spc_tmp = [1, 1, 1];
img_CT = double(img_CT);
img_CT = img_thr(img_CT, -1024, 0, 1);
img_GT_cur = img_GT(:,:,:,1);
img_GT_cur = imresize3(img_GT_cur, round(bszv .* units .* spc_tmp), 'linear');
GT_size = size(img_GT_cur);
img_CT_rsz = volresize(img_CT, GT_size, 1);
img_CT_s = squeeze(img_CT_rsz(:,:,index_slice));
img_CT_s = imrotate(img_CT_s,90);
CT_RGB = repmat(img_CT_s,[1,1,3]);
img_GT = double(img_GT);
mask = double(mask);

% img_GT = volresize(img_GT, round(bszv .* units .* spc_tmp), 1);
img_GT_s_mask_all = [];
for num_fig = 1:num_phase
    img_GT_cur = img_GT(:,:,:,num_fig);
    img_GT_cur = imresize3(img_GT_cur, round(bszv .* units .* spc_tmp), 'linear');
    mask = imresize3(mask, round(bszv .* units .* spc_tmp), 'nearest');
    img_GT_mask = img_GT_cur.*mask;

    new_sz = size(img_GT_cur);

    
    
    
    img_GT_s = squeeze(img_GT_mask(:,:,index_slice));
    img_GT_s = imrotate(img_GT_s,90);
    
    mask_s = squeeze(double(mask(:,:,index_slice)));
    mask_s = imrotate(mask_s,90);
    img_GT_s_mask_cur = double(img_GT_s).*double(mask_s);
    img_GT_s_mask_all = cat(3,img_GT_s_mask_all,img_GT_s_mask_cur);
end

CT_alpha = 1 - mask_s;


max_value_slice = max(img_GT_s_mask_all(:));
min_value_slice = min(img_GT_s_mask_all(:));
img_GT_s_mask_all = img_thr(img_GT_s_mask_all, min_value_slice, max_value_slice, 1);
u_value = unique(img_GT_s_mask_all);
dx = 1/length(u_value);
inter = [0:dx:1];
for un_idx = 1:length(u_value)
    img_GT_s_mask_all(find(img_GT_s_mask_all==u_value(un_idx))) = inter(un_idx);
end

for num_fig = 1:num_phase
    figure
    img_GT_s_mask_cur = img_GT_s_mask_all(:,:,num_fig);
    if strcmp(colormap_name,'jet')
        imshow(img_GT_s_mask_cur,'Colormap',jet(255));
    else
        if strcmp(colormap_name,'parula')
            imshow(img_GT_s_mask_cur,'Colormap',parula(255));
        else
            if strcmp(colormap_name,'hot')
                imshow(img_GT_s_mask_cur,'Colormap',hot(1000));
            else
                if strcmp(colormap_name,'autumn')
                    imshow(img_GT_s_mask_cur,'Colormap',autumn(size(u_value,1)));
                end
            end
        end
    end
    hold on;
    h = imshow(CT_RGB,[]);
    set(h, 'alphadata', CT_alpha);
end
end
