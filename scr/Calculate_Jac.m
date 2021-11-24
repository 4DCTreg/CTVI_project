function [Tptv_rsz_new] = Calculate_Jac(vol_ex, Tptv_rsz, spc)
im_sz = size(vol_ex);

Tptv_rsz_new = Tptv_rsz;

Tx = padarray(Tptv_rsz(:,:,:,1),[1, 1, 1]);
Ty = padarray(Tptv_rsz(:,:,:,2),[1, 1, 1]);
Tz = padarray(Tptv_rsz(:,:,:,3),[1, 1, 1]);
inv_index_matrix = zeros(im_sz);
top_co = 1;

parfor index = 1:im_sz(1)*im_sz(2)*im_sz(3)
    
    [ix, iy, iz] = ind2sub(im_sz, index);
    ix_add = ix+1;
    iy_add = iy+1;
    iz_add = iz+1;
    
    tporpx = Tx(ix_add,iy_add,iz_add);
    tporpy = Ty(ix_add,iy_add,iz_add);
    tporpz = Tz(ix_add,iy_add,iz_add);
    
    or_p_up_x = Tx(ix_add,iy_add,iz_add-1);
    or_p_up_y = Ty(ix_add,iy_add,iz_add-1);
    or_p_up_z = Tz(ix_add,iy_add,iz_add-1) + spc(3);
    
    or_p_down_x = Tx(ix_add,iy_add,iz_add+1);
    or_p_down_y = Ty(ix_add,iy_add,iz_add+1);
    or_p_down_z = Tz(ix_add,iy_add,iz_add+1) - spc(3);
    
    or_p_left_x = Tx(ix_add,iy_add-1,iz_add);
    or_p_left_y = Ty(ix_add,iy_add-1,iz_add) - spc(2);
    or_p_left_z = Tz(ix_add,iy_add-1,iz_add);
    
    or_p_right_x = Tx(ix_add,iy_add+1,iz_add);
    or_p_right_y = Ty(ix_add,iy_add+1,iz_add) + spc(2);
    or_p_right_z = Tz(ix_add,iy_add+1,iz_add);
    
    or_p_anter_x = Tx(ix_add+1,iy_add,iz_add) - spc(1);
    or_p_anter_y = Ty(ix_add+1,iy_add,iz_add);
    or_p_anter_z = Tz(ix_add+1,iy_add,iz_add);
    
    or_p_poster_x = Tx(ix_add-1,iy_add,iz_add) + spc(1);
    or_p_poster_y = Ty(ix_add-1,iy_add,iz_add);
    or_p_poster_z = Tz(ix_add-1,iy_add,iz_add);
    
    %Jaco1 point 1 4 5
    %Jaco_1=(p_right_x-orpx)*(orpy-p_anter_y)*(p_up_z-orpz)+(orpx-p_anter_x)*(p_up_y-orpy)*(p_right_z-orpz)+(p_right_y-orpy)*(orpz-p_anter_z)*(p_up_x-orpx)-(p_right_z-orpz)*(orpy-p_anter_y)*(p_up_x-orpx)-(p_right_y-orpy)*(orpx-p_anter_x)*(p_up_z-orpz)-(orpz-p_anter_z)*(p_up_y-orpy)*(p_right_x-orpx);
    Jaco_1 = (tporpx-or_p_anter_x)*(or_p_right_y-tporpy)*(or_p_up_z-tporpz)+(tporpy-or_p_anter_y)*(or_p_right_z-tporpz)*(or_p_up_x-tporpx)+(or_p_right_x-tporpx)*(or_p_up_y-tporpy)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(or_p_right_y-tporpy)*(or_p_up_x-tporpx)-(tporpy-or_p_anter_y)*(or_p_right_x-tporpx)*(or_p_up_z-tporpz)-(or_p_right_z-tporpz)*(or_p_up_y-tporpy)*(tporpx-or_p_anter_x);
    Jaco_1 = Jaco_1*top_co;
    
    
    
    %Jaco2 point 1 3 5
    %Jaco_2=(orpx-p_left_x)*(orpy-p_anter_y)*(p_up_z-orpz)+(orpx-p_anter_x)*(p_up_y-orpy)*(orpz-p_left_z)+(orpy-p_left_y)*(orpz-p_anter_z)*(p_up_x-orpx)-(orpz-p_left_z)*(orpy-p_anter_y)*(p_up_x-orpx)-(orpy-p_left_y)*(orpx-p_anter_x)*(p_up_z-orpz)-(orpz-p_anter_z)*(p_up_y-orpy)*(orpx-p_left_x);
    Jaco_2 = (tporpx-or_p_anter_x)*(tporpy-or_p_left_y)*(or_p_up_z-tporpz)+(tporpy-or_p_anter_y)*(tporpz-or_p_left_z)*(or_p_up_x-tporpx)+(tporpx-or_p_left_x)*(or_p_up_y-tporpy)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(tporpy-or_p_left_y)*(or_p_up_x-tporpx)-(tporpy-or_p_anter_y)*(tporpx-or_p_left_x)*(or_p_up_z-tporpz)-(tporpz-or_p_left_z)*(or_p_up_y-tporpy)*(tporpx-or_p_anter_x);
    Jaco_2 = Jaco_2*top_co;
    
    
    %Jaco3 point 2 4 5
    %Jaco_3=(p_right_x-orpx)*(orpy-p_anter_y)*(orpz-p_down_z)+(orpx-p_anter_x)*(orpy-p_down_y)*(p_right_z-orpz)+(p_right_y-orpy)*(orpz-p_anter_z)*(orpx-p_down_x)-(p_right_z-orpz)*(orpy-p_anter_y)*(orpx-p_down_x)-(p_right_y-orpy)*(orpx-p_anter_x)*(orpz-p_down_z)-(orpz-p_anter_z)*(orpy-p_down_y)*(p_right_x-orpx);
    Jaco_3 = (tporpx-or_p_anter_x)*(or_p_right_y-tporpy)*(tporpz-or_p_down_z)+(tporpy-or_p_anter_y)*(or_p_right_z-tporpz)*(tporpx-or_p_down_x)+(or_p_right_x-tporpx)*(tporpy-or_p_down_y)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(or_p_right_y-tporpy)*(tporpx-or_p_down_x)-(tporpy-or_p_anter_y)*(or_p_right_x-tporpx)*(tporpz-or_p_down_z)-(or_p_right_z-tporpz)*(tporpy-or_p_down_y)*(tporpx-or_p_anter_x);
    Jaco_3 = Jaco_3*top_co;
    
    
    %Jaco4 point 2 3 5
    %Jaco_4=(orpx-p_left_x)*(orpy-p_anter_y)*(orpz-p_down_z)+(orpx-p_anter_x)*(orpy-p_down_y)*(orpz-p_left_z)+(orpy-p_left_y)*(orpz-p_anter_z)*(orpx-p_down_x)-(orpz-p_left_z)*(orpy-p_anter_y)*(orpx-p_down_x)-(orpy-p_left_y)*(orpx-p_anter_x)*(orpz-p_down_z)-(orpz-p_anter_z)*(orpy-p_down_y)*(orpx-p_left_x);
    Jaco_4 = (tporpx-or_p_anter_x)*(tporpy-or_p_left_y)*(tporpz-or_p_down_z)+(tporpy-or_p_anter_y)*(tporpz-or_p_left_z)*(tporpx-or_p_down_x)+(tporpx-or_p_left_x)*(tporpy-or_p_down_y)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(tporpy-or_p_left_y)*(tporpx-or_p_down_x)-(tporpy-or_p_anter_y)*(tporpx-or_p_left_x)*(tporpz-or_p_down_z)-(tporpz-or_p_left_z)*(tporpy-or_p_down_y)*(tporpx-or_p_anter_x);
    Jaco_4 = Jaco_4*top_co;
    
    
    %Jaco5 point 1 4 6 fff
    %Jaco_5=(p_right_x-orpx)*(p_poster_y-orpy)*(p_up_z-orpz)+(p_poster_x-orpx)*(p_up_y-orpy)*(p_right_z-orpz)+(p_right_y-orpy)*(p_poster_z-orpz)*(p_up_x-orpx)-(p_right_z-orpz)*(p_poster_y-orpy)*(p_up_x-orpx)-(p_right_y-orpy)*(p_poster_x-orpx)*(p_up_z-orpz)-(p_poster_z-orpz)*(p_up_y-orpy)*(p_right_x-orpx);
    Jaco_5 = (or_p_poster_x-tporpx)*(or_p_right_y-tporpy)*(or_p_up_z-tporpz)+(or_p_poster_y-tporpy)*(or_p_right_z-tporpz)*(or_p_up_x-tporpx)+(or_p_right_x-tporpx)*(or_p_up_y-tporpy)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(or_p_right_y-tporpy)*(or_p_up_x-tporpx)-(or_p_poster_y-tporpy)*(or_p_right_x-tporpx)*(or_p_up_z-tporpz)-(or_p_right_z-tporpz)*(or_p_up_y-tporpy)*(or_p_poster_x-tporpx);
    Jaco_5 = Jaco_5*top_co;
    
    
    
    %Jaco6 point 2 3 6
    %Jaco_6=(orpx-p_left_x)*(p_poster_y-orpy)*(orpz-p_down_z)+(p_poster_x-orpx)*(orpy-p_down_y)*(orpz-p_left_z)+(orpy-p_left_y)*(p_poster_z-orpz)*(orpx-p_down_x)-(orpz-p_left_z)*(p_poster_y-orpy)*(orpx-p_down_x)-(orpy-p_left_y)*(p_poster_x-orpx)*(orpz-p_down_z)-(p_poster_z-orpz)*(orpy-p_down_y)*(orpx-p_left_x);
    Jaco_6 = (or_p_poster_x-tporpx)*(tporpy-or_p_left_y)*(tporpz-or_p_down_z)+(or_p_poster_y-tporpy)*(tporpz-or_p_left_z)*(tporpx-or_p_down_x)+(tporpx-or_p_left_x)*(tporpy-or_p_down_y)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(tporpy-or_p_left_y)*(tporpx-or_p_down_x)-(or_p_poster_y-tporpy)*(tporpx-or_p_left_x)*(tporpz-or_p_down_z)-(tporpz-or_p_left_z)*(tporpy-or_p_down_y)*(or_p_poster_x-tporpx);
    Jaco_6 = Jaco_6*top_co;
    
    
    %Jaco7 point 2 4 6
    %Jaco_7=(p_right_x-orpx)*(p_poster_y-orpy)*(orpz-p_down_z)+(p_poster_x-orpx)*(orpy-p_down_y)*(p_right_z-orpz)+(p_right_y-orpy)*(p_poster_z-orpz)*(orpx-p_down_x)-(p_right_z-orpz)*(p_poster_y-orpy)*(orpx-p_down_x)-(p_right_y-orpy)*(p_poster_x-orpx)*(orpz-p_down_z)-(p_poster_z-orpz)*(orpy-p_down_y)*(-p_right_x-orpx);
    
    Jaco_7 = (or_p_poster_x-tporpx)*(or_p_right_y-tporpy)*(tporpz-or_p_down_z)+(or_p_poster_y-tporpy)*(or_p_right_z-tporpz)*(tporpx-or_p_down_x)+(or_p_right_x-tporpx)*(tporpy-or_p_down_y)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(or_p_right_y-tporpy)*(tporpx-or_p_down_x)-(or_p_poster_y-tporpy)*(or_p_right_x-tporpx)*(tporpz-or_p_down_z)-(or_p_right_z-tporpz)*(tporpy-or_p_down_y)*(or_p_poster_x-tporpx);
    Jaco_7 = Jaco_7*top_co;
    
    
    %Jaco8 point 1 3 6
    %Jaco_8=(orpx-p_left_x)*(p_poster_y-orpy)*(p_up_z-orpz)+(p_poster_x-orpx)*(p_up_y-orpy)*(orpz-p_left_z)+(orpy-p_left_y)*(p_poster_z-orpz)*(p_up_x-orpx)-(orpz-p_left_z)*(p_poster_y-orpy)*(p_up_x-orpx)-(orpy-p_left_y)*(p_poster_x-orpx)*(p_up_z-orpz)-(p_poster_z-orpz)*(p_up_y-orpy)*(orpx-p_left_x);
    Jaco_8 = (or_p_poster_x-tporpx)*(tporpy-or_p_left_y)*(or_p_up_z-tporpz)+(or_p_poster_y-tporpy)*(tporpz-or_p_left_z)*(or_p_up_x-tporpx)+(tporpx-or_p_left_x)*(or_p_up_y-tporpy)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(tporpy-or_p_left_y)*(or_p_up_x-tporpx)-(or_p_poster_y-tporpy)*(tporpx-or_p_left_x)*(or_p_up_z-tporpz)-(tporpz-or_p_left_z)*(or_p_up_y-tporpy)*(or_p_poster_x-tporpx);
    Jaco_8 = Jaco_8*top_co;
    
    %     vent_img(index) = mean([Jaco_1,Jaco_2,Jaco_3,Jaco_4,Jaco_5,Jaco_6,Jaco_7,Jaco_8]);
    current_jac = [Jaco_1,Jaco_2,Jaco_3,Jaco_4,Jaco_5,Jaco_6,Jaco_7,Jaco_8];
    inv_index = ~isempty(find(current_jac<0));
    inv_index_matrix(index) = inv_index;
    
    
end


Tx_smooth = imgaussfilt3(Tptv_rsz(:,:,:,1),1,'FilterSize',9);
Ty_smooth = imgaussfilt3(Tptv_rsz(:,:,:,2),1,'FilterSize',9);
Tz_smooth = imgaussfilt3(Tptv_rsz(:,:,:,3),1,'FilterSize',9);

sm = find(inv_index_matrix==1);
Tx_or = Tptv_rsz(:,:,:,1);
Ty_or = Tptv_rsz(:,:,:,2);
Tz_or = Tptv_rsz(:,:,:,3);
Tx_or(sm) = Tx_smooth(sm);
Ty_or(sm) = Ty_smooth(sm);
Tz_or(sm) = Tz_smooth(sm);
Tptv_rsz_new(:,:,:,1) = Tx_or;
Tptv_rsz_new(:,:,:,2) = Ty_or;
Tptv_rsz_new(:,:,:,3) = Tz_or;
end