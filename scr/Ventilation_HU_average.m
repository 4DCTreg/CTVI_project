function vent_img = Ventilation_HU_average(vol_ct)
vol_ct = double(vol_ct);
[r,c,s,N] = size(vol_ct);
vent_img = zeros(r,c,s);
for id = 1:N
    CT_cur = vol_ct(:,:,:,id);
    vent_img_cur = (-CT_cur/1000).*((CT_cur+1000)/1000);
    vent_img = vent_img+vent_img_cur;
     
end
vent_img = vent_img./6;