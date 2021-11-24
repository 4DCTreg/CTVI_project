function crop_v = crop_mask(avg_mask,info,d)
sz_ct = info;

for r_index = 1:sz_ct(1)
    slice_cur = squeeze(avg_mask(r_index,:,:));
    min_crop_r = r_index;
    if any(slice_cur(:))
        break
    end
end

for r_index = sz_ct(1):-1:1
    slice_cur = squeeze(avg_mask(r_index,:,:));
    if any(slice_cur(:))
        max_crop_r = r_index;
        break
    end
end

for c_index = 1:sz_ct(2)
    slice_cur = squeeze(avg_mask(:,c_index,:));
    if any(slice_cur(:))
        min_crop_c = c_index;
        break
    end
end

for c_index = sz_ct(2):-1:1
    slice_cur = squeeze(avg_mask(:,c_index,:));
    if any(slice_cur(:))
        max_crop_c = c_index;
        break
    end
end

for s_index = 1:sz_ct(3)
    slice_cur = squeeze(avg_mask(:,:,s_index));
    if any(slice_cur(:))
        min_crop_s = s_index;
        break
    end
end

for s_index = sz_ct(3):-1:1
    slice_cur = squeeze(avg_mask(:,:,s_index));
    if any(slice_cur(:))
        max_crop_s = s_index;
        break
    end
end

crop_v = [ max(1, min_crop_r - d(1)), min(sz_ct(1), max_crop_r + d(1)); ... 
           max(1, min_crop_c - d(2)), min(sz_ct(2), max_crop_c + d(2)); ... 
           max(1, min_crop_s - d(3)), min(sz_ct(3), max_crop_s + d(3));];
end