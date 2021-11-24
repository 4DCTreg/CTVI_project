function [new_crop,min_h_distance] =mulit_sup_corp(crop_v,sup_real_sz,h_distance)
new_crop = zeros(3,2);
min_h_distance = zeros(3,2);
if crop_v(1,1)-h_distance>0
    new_crop(1,1) = crop_v(1,1)-h_distance;
    min_h_distance(1,1) = h_distance;
else
    new_crop(1,1) = 1;
    min_h_distance(1,1) = crop_v(1,1)-1;
end

if sup_real_sz(1)-h_distance>=crop_v(1,2)
    new_crop(1,2) = crop_v(1,2)+h_distance;
    min_h_distance(1,2) = h_distance;
else
    new_crop(1,2) = sup_real_sz(1);
    min_h_distance(1,2) = sup_real_sz(1)-crop_v(1,2);
end



if crop_v(2,1)-h_distance>0
    new_crop(2,1) = crop_v(2,1)-h_distance;
    min_h_distance(2,1) = h_distance;
else
    new_crop(2,1) = 1;
    min_h_distance(2,1) = crop_v(2,1)-1;
end

if sup_real_sz(2)-h_distance>=crop_v(2,2)
    new_crop(2,2) = crop_v(2,2)+h_distance;
    min_h_distance(2,2) = h_distance;
else
    new_crop(2,2) = sup_real_sz(2);
    min_h_distance(2,2) = sup_real_sz(2)-crop_v(2,2);
end



if crop_v(3,1)-h_distance>0
    new_crop(3,1) = crop_v(3,1)-h_distance;
    min_h_distance(3,1) = h_distance;
else
    new_crop(3,1) = 1;
    min_h_distance(3,1) = crop_v(3,1)-1;
end

if sup_real_sz(3)-h_distance>=crop_v(3,2)
    new_crop(3,2) = crop_v(3,2)+h_distance;
    min_h_distance(3,2) = h_distance;
else
    new_crop(3,2) = sup_real_sz(3);
    min_h_distance(3,2) = sup_real_sz(3)-crop_v(3,2);
end

min_h_distance = min(min_h_distance,[],2);
end
    