clear;
load('HOMRF_jac_debbug.mat')
im_sz = size(vol_ex);
% vent_img = zeros(im_sz);
% vent_jac = cell(im_sz);
vent_jac = zeros(im_sz(1)*im_sz(2)*im_sz(3),8);
% vent_jac = zeros(im_sz(1)*im_sz(2)*im_sz(3),8);

Tx = padarray(Tptv_rsz(:,:,:,1),[1, 1, 1]);
Ty = padarray(Tptv_rsz(:,:,:,2),[1, 1, 1]);
Tz = padarray(Tptv_rsz(:,:,:,3),[1, 1, 1]);

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
    if Jaco_1 <= 0
        Jaco_1 = 1;
    end
    
    
    
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
%     vent_img(index) = mean(log(Jaco_1));
    vent_jac(index,:) = [Jaco_1,Jaco_2,Jaco_3,Jaco_4,Jaco_5,Jaco_6,Jaco_7,Jaco_8];
    
    
end
% vent_img(find(vent_img)<0) = 0;
% 3D linear interpolation
% divide into 8 groups, each groups has 8 points
% groups
%   1！4
%  /  /
% 2！3
%   5！8
%  /  /
% 6！7
if mod(im_sz(1),2)==0
    g_1_1_x = [1:2:im_sz(1)];
    g_1_2_x = g_1_1_x+1;
    g_1_3_x = g_1_2_x;
    g_1_4_x = g_1_1_x;
    g_1_5_x = g_1_1_x;
    g_1_6_x = g_1_2_x;
    g_1_7_x = g_1_2_x;
    g_1_8_x = g_1_1_x;
    
%     g_2_1_x = g_1_1_x+1;
    g_2_1_x = [2:2:im_sz(1)-1];
    g_2_2_x = g_2_1_x+1;
    g_2_3_x = g_2_2_x;
    g_2_4_x = g_2_1_x;
    g_2_5_x = g_2_1_x;
    g_2_6_x = g_2_2_x;
    g_2_7_x = g_2_2_x;
    g_2_8_x = g_2_1_x;
    
    g_3_1_x = g_2_1_x;
    g_3_2_x = g_3_1_x+1;
    g_3_3_x = g_3_2_x;
    g_3_4_x = g_3_1_x;
    g_3_5_x = g_3_1_x;
    g_3_6_x = g_3_2_x;
    g_3_7_x = g_3_2_x;
    g_3_8_x = g_3_1_x;
    
    g_4_1_x = g_1_1_x;
    g_4_2_x = g_4_1_x+1;
    g_4_3_x = g_4_2_x;
    g_4_4_x = g_4_1_x;
    g_4_5_x = g_4_1_x;
    g_4_6_x = g_4_2_x;
    g_4_7_x = g_4_2_x;
    g_4_8_x = g_4_1_x;
    
    g_5_1_x = g_1_1_x;
    g_5_2_x = g_5_1_x+1;
    g_5_3_x = g_5_2_x;
    g_5_4_x = g_5_1_x;
    g_5_5_x = g_5_1_x;
    g_5_6_x = g_5_2_x;
    g_5_7_x = g_5_2_x;
    g_5_8_x = g_5_1_x;
    
    g_6_1_x = g_2_1_x;
    g_6_2_x = g_6_1_x+1;
    g_6_3_x = g_6_2_x;
    g_6_4_x = g_6_1_x;
    g_6_5_x = g_6_1_x;
    g_6_6_x = g_6_2_x;
    g_6_7_x = g_6_2_x;
    g_6_8_x = g_6_1_x;
    
    g_7_1_x = g_2_1_x;
    g_7_2_x = g_7_1_x+1;
    g_7_3_x = g_7_2_x;
    g_7_4_x = g_7_1_x;
    g_7_5_x = g_7_1_x;
    g_7_6_x = g_7_2_x;
    g_7_7_x = g_7_2_x;
    g_7_8_x = g_7_1_x;
    
    g_8_1_x = g_1_1_x;
    g_8_2_x = g_8_1_x+1;
    g_8_3_x = g_8_2_x;
    g_8_4_x = g_8_1_x;
    g_8_5_x = g_8_1_x;
    g_8_6_x = g_8_2_x;
    g_8_7_x = g_8_2_x;
    g_8_8_x = g_8_1_x;
else
    g_1_1_x = [1:2:im_sz(1)-1];
    g_1_2_x = g_1_1_x+1;
    g_1_3_x = g_1_2_x;
    g_1_4_x = g_1_1_x;
    g_1_5_x = g_1_1_x;
    g_1_6_x = g_1_2_x;
    g_1_7_x = g_1_2_x;
    g_1_8_x = g_1_1_x;
    
%     g_2_1_x = g_1_1_x+1;
    g_2_1_x = [2:2:im_sz(1)-1];
    g_2_2_x = g_2_1_x+1;
    g_2_3_x = g_2_2_x;
    g_2_4_x = g_2_1_x;
    g_2_5_x = g_2_1_x;
    g_2_6_x = g_2_2_x;
    g_2_7_x = g_2_2_x;
    g_2_8_x = g_2_1_x;
    
    g_3_1_x = g_2_1_x;
    g_3_2_x = g_3_1_x+1;
    g_3_3_x = g_3_2_x;
    g_3_4_x = g_3_1_x;
    g_3_5_x = g_3_1_x;
    g_3_6_x = g_3_2_x;
    g_3_7_x = g_3_2_x;
    g_3_8_x = g_3_1_x;
    
    g_4_1_x = g_1_1_x;
    g_4_2_x = g_4_1_x+1;
    g_4_3_x = g_4_2_x;
    g_4_4_x = g_4_1_x;
    g_4_5_x = g_4_1_x;
    g_4_6_x = g_4_2_x;
    g_4_7_x = g_4_2_x;
    g_4_8_x = g_4_1_x;
    
    g_5_1_x = g_1_1_x;
    g_5_2_x = g_5_1_x+1;
    g_5_3_x = g_5_2_x;
    g_5_4_x = g_5_1_x;
    g_5_5_x = g_5_1_x;
    g_5_6_x = g_5_2_x;
    g_5_7_x = g_5_2_x;
    g_5_8_x = g_5_1_x;
    
    g_6_1_x = g_2_1_x;
    g_6_2_x = g_6_1_x+1;
    g_6_3_x = g_6_2_x;
    g_6_4_x = g_6_1_x;
    g_6_5_x = g_6_1_x;
    g_6_6_x = g_6_2_x;
    g_6_7_x = g_6_2_x;
    g_6_8_x = g_6_1_x;
    
    g_7_1_x = g_2_1_x;
    g_7_2_x = g_7_1_x+1;
    g_7_3_x = g_7_2_x;
    g_7_4_x = g_7_1_x;
    g_7_5_x = g_7_1_x;
    g_7_6_x = g_7_2_x;
    g_7_7_x = g_7_2_x;
    g_7_8_x = g_7_1_x;
    
    g_8_1_x = g_1_1_x;
    g_8_2_x = g_8_1_x+1;
    g_8_3_x = g_8_2_x;
    g_8_4_x = g_8_1_x;
    g_8_5_x = g_8_1_x;
    g_8_6_x = g_8_2_x;
    g_8_7_x = g_8_2_x;
    g_8_8_x = g_8_1_x;
end

if mod(im_sz(2),2) == 0
    g_1_1_y = [1:2:im_sz(2)];
    g_1_2_y = g_1_1_y;
    g_1_3_y = g_1_1_y+1;
    g_1_4_y = g_1_3_y;
    g_1_5_y = g_1_1_y;
    g_1_6_y = g_1_1_y;
    g_1_7_y = g_1_3_y;
    g_1_8_y = g_1_4_y;
    
    g_2_1_y = g_1_1_y;
    g_2_2_y = g_2_1_y;
    g_2_3_y = g_2_1_y+1;
    g_2_4_y = g_2_3_y;
    g_2_5_y = g_2_1_y;
    g_2_6_y = g_2_1_y;
    g_2_7_y = g_2_3_y;
    g_2_8_y = g_2_4_y;
    
%     g_3_1_y = g_1_1_y+1;
    g_3_1_y = [2:2:im_sz(2)-1];
    g_3_2_y = g_3_1_y;
    g_3_3_y = g_3_1_y+1;
    g_3_4_y = g_3_3_y;
    g_3_5_y = g_3_1_y;
    g_3_6_y = g_3_1_y;
    g_3_7_y = g_3_3_y;
    g_3_8_y = g_3_4_y;
    
    g_4_1_y = g_3_1_y;
    g_4_2_y = g_4_1_y;
    g_4_3_y = g_4_1_y+1;
    g_4_4_y = g_4_3_y;
    g_4_5_y = g_4_1_y;
    g_4_6_y = g_4_1_y;
    g_4_7_y = g_4_3_y;
    g_4_8_y = g_4_4_y;
    
    g_5_1_y = g_1_1_y;
    g_5_2_y = g_5_1_y;
    g_5_3_y = g_5_1_y+1;
    g_5_4_y = g_5_3_y;
    g_5_5_y = g_5_1_y;
    g_5_6_y = g_5_1_y;
    g_5_7_y = g_5_3_y;
    g_5_8_y = g_5_4_y;
    
    g_6_1_y = g_1_1_y;
    g_6_2_y = g_6_1_y;
    g_6_3_y = g_6_1_y+1;
    g_6_4_y = g_6_3_y;
    g_6_5_y = g_6_1_y;
    g_6_6_y = g_6_1_y;
    g_6_7_y = g_6_3_y;
    g_6_8_y = g_6_4_y;
    
    g_7_1_y = g_3_1_y;
    g_7_2_y = g_7_1_y;
    g_7_3_y = g_7_1_y+1;
    g_7_4_y = g_7_3_y;
    g_7_5_y = g_7_1_y;
    g_7_6_y = g_7_1_y;
    g_7_7_y = g_7_3_y;
    g_7_8_y = g_7_4_y;
    
    g_8_1_y = g_3_1_y;
    g_8_2_y = g_8_1_y;
    g_8_3_y = g_8_1_y+1;
    g_8_4_y = g_8_3_y;
    g_8_5_y = g_8_1_y;
    g_8_6_y = g_8_1_y;
    g_8_7_y = g_8_3_y;
    g_8_8_y = g_8_4_y;
    
else
    g_1_1_y = [1:2:im_sz(2)-1];
    g_1_2_y = g_1_1_y;
    g_1_3_y = g_1_1_y+1;
    g_1_4_y = g_1_3_y;
    g_1_5_y = g_1_1_y;
    g_1_6_y = g_1_1_y;
    g_1_7_y = g_1_3_y;
    g_1_8_y = g_1_4_y;
    
    g_2_1_y = g_1_1_y;
    g_2_2_y = g_2_1_y;
    g_2_3_y = g_2_1_y+1;
    g_2_4_y = g_2_3_y;
    g_2_5_y = g_2_1_y;
    g_2_6_y = g_2_1_y;
    g_2_7_y = g_2_3_y;
    g_2_8_y = g_2_4_y;
    
    g_3_1_y = [2:2:im_sz(2)-1];
    g_3_2_y = g_3_1_y;
    g_3_3_y = g_3_1_y+1;
    g_3_4_y = g_3_3_y;
    g_3_5_y = g_3_1_y;
    g_3_6_y = g_3_1_y;
    g_3_7_y = g_3_3_y;
    g_3_8_y = g_3_4_y;
    
    g_4_1_y = g_3_1_y;
    g_4_2_y = g_4_1_y;
    g_4_3_y = g_4_1_y+1;
    g_4_4_y = g_4_3_y;
    g_4_5_y = g_4_1_y;
    g_4_6_y = g_4_1_y;
    g_4_7_y = g_4_3_y;
    g_4_8_y = g_4_4_y;
    
    g_5_1_y = g_1_1_y;
    g_5_2_y = g_5_1_y;
    g_5_3_y = g_5_1_y+1;
    g_5_4_y = g_5_3_y;
    g_5_5_y = g_5_1_y;
    g_5_6_y = g_5_1_y;
    g_5_7_y = g_5_3_y;
    g_5_8_y = g_5_4_y;
    
    g_6_1_y = g_1_1_y;
    g_6_2_y = g_6_1_y;
    g_6_3_y = g_6_1_y+1;
    g_6_4_y = g_6_3_y;
    g_6_5_y = g_6_1_y;
    g_6_6_y = g_6_1_y;
    g_6_7_y = g_6_3_y;
    g_6_8_y = g_6_4_y;
    
    g_7_1_y = g_3_1_y;
    g_7_2_y = g_7_1_y;
    g_7_3_y = g_7_1_y+1;
    g_7_4_y = g_7_3_y;
    g_7_5_y = g_7_1_y;
    g_7_6_y = g_7_1_y;
    g_7_7_y = g_7_3_y;
    g_7_8_y = g_7_4_y;
    
    g_8_1_y = g_3_1_y;
    g_8_2_y = g_8_1_y;
    g_8_3_y = g_8_1_y+1;
    g_8_4_y = g_8_3_y;
    g_8_5_y = g_8_1_y;
    g_8_6_y = g_8_1_y;
    g_8_7_y = g_8_3_y;
    g_8_8_y = g_8_4_y;
end

if mod(im_sz(3),2) == 0
    g_1_1_z = [1:2:im_sz(3)];
    g_1_2_z = g_1_1_z;
    g_1_3_z = g_1_1_z;
    g_1_4_z = g_1_1_z;
    g_1_5_z = g_1_1_z+1;
    g_1_6_z = g_1_5_z;
    g_1_7_z = g_1_5_z;
    g_1_8_z = g_1_5_z;
    
    g_2_1_z = g_1_1_z;
    g_2_2_z = g_2_1_z;
    g_2_3_z = g_2_1_z;
    g_2_4_z = g_2_1_z;
    g_2_5_z = g_2_1_z+1;
    g_2_6_z = g_2_5_z;
    g_2_7_z = g_2_5_z;
    g_2_8_z = g_2_5_z;
    
    g_3_1_z = g_1_1_z;
    g_3_2_z = g_3_1_z;
    g_3_3_z = g_3_1_z;
    g_3_4_z = g_3_1_z;
    g_3_5_z = g_3_1_z+1;
    g_3_6_z = g_3_5_z;
    g_3_7_z = g_3_5_z;
    g_3_8_z = g_3_5_z;
    
    g_4_1_z = g_1_1_z;
    g_4_2_z = g_4_1_z;
    g_4_3_z = g_4_1_z;
    g_4_4_z = g_4_1_z;
    g_4_5_z = g_4_1_z+1;
    g_4_6_z = g_4_5_z;
    g_4_7_z = g_4_5_z;
    g_4_8_z = g_4_5_z;
    
    g_5_1_z = [2:2:im_sz(3)-1];
    g_5_2_z = g_5_1_z;
    g_5_3_z = g_5_1_z;
    g_5_4_z = g_5_1_z;
    g_5_5_z = g_5_1_z+1;
    g_5_6_z = g_5_5_z;
    g_5_7_z = g_5_5_z;
    g_5_8_z = g_5_5_z;
    
    g_6_1_z = g_5_1_z;
    g_6_2_z = g_6_1_z;
    g_6_3_z = g_6_1_z;
    g_6_4_z = g_6_1_z;
    g_6_5_z = g_6_1_z+1;
    g_6_6_z = g_6_5_z;
    g_6_7_z = g_6_5_z;
    g_6_8_z = g_6_5_z;
    
    g_7_1_z = g_5_1_z;
    g_7_2_z = g_7_1_z;
    g_7_3_z = g_7_1_z;
    g_7_4_z = g_7_1_z;
    g_7_5_z = g_7_1_z+1;
    g_7_6_z = g_7_5_z;
    g_7_7_z = g_7_5_z;
    g_7_8_z = g_7_5_z;
    
    g_8_1_z = g_5_1_z;
    g_8_2_z = g_8_1_z;
    g_8_3_z = g_8_1_z;
    g_8_4_z = g_8_1_z;
    g_8_5_z = g_8_1_z+1;
    g_8_6_z = g_8_5_z;
    g_8_7_z = g_8_5_z;
    g_8_8_z = g_8_5_z;
else
    g_1_1_z = [1:2:im_sz(3)-1];
    g_1_2_z = g_1_1_z;
    g_1_3_z = g_1_1_z;
    g_1_4_z = g_1_1_z;
    g_1_5_z = g_1_1_z+1;
    g_1_6_z = g_1_5_z;
    g_1_7_z = g_1_5_z;
    g_1_8_z = g_1_5_z;
    
    g_2_1_z = g_1_1_z;
    g_2_2_z = g_2_1_z;
    g_2_3_z = g_2_1_z;
    g_2_4_z = g_2_1_z;
    g_2_5_z = g_2_1_z+1;
    g_2_6_z = g_2_5_z;
    g_2_7_z = g_2_5_z;
    g_2_8_z = g_2_5_z;
    
    g_3_1_z = g_1_1_z;
    g_3_2_z = g_3_1_z;
    g_3_3_z = g_3_1_z;
    g_3_4_z = g_3_1_z;
    g_3_5_z = g_3_1_z+1;
    g_3_6_z = g_3_5_z;
    g_3_7_z = g_3_5_z;
    g_3_8_z = g_3_5_z;
    
    g_4_1_z = g_1_1_z;
    g_4_2_z = g_4_1_z;
    g_4_3_z = g_4_1_z;
    g_4_4_z = g_4_1_z;
    g_4_5_z = g_4_1_z+1;
    g_4_6_z = g_4_5_z;
    g_4_7_z = g_4_5_z;
    g_4_8_z = g_4_5_z;
    
    g_5_1_z = [2:2:im_sz(3)-1];
    g_5_2_z = g_5_1_z;
    g_5_3_z = g_5_1_z;
    g_5_4_z = g_5_1_z;
    g_5_5_z = g_5_1_z+1;
    g_5_6_z = g_5_5_z;
    g_5_7_z = g_5_5_z;
    g_5_8_z = g_5_5_z;
    
    g_6_1_z = g_5_1_z;
    g_6_2_z = g_6_1_z;
    g_6_3_z = g_6_1_z;
    g_6_4_z = g_6_1_z;
    g_6_5_z = g_6_1_z+1;
    g_6_6_z = g_6_5_z;
    g_6_7_z = g_6_5_z;
    g_6_8_z = g_6_5_z;
    
    g_7_1_z = g_5_1_z;
    g_7_2_z = g_7_1_z;
    g_7_3_z = g_7_1_z;
    g_7_4_z = g_7_1_z;
    g_7_5_z = g_7_1_z+1;
    g_7_6_z = g_7_5_z;
    g_7_7_z = g_7_5_z;
    g_7_8_z = g_7_5_z;
    
    g_8_1_z = g_5_1_z;
    g_8_2_z = g_8_1_z;
    g_8_3_z = g_8_1_z;
    g_8_4_z = g_8_1_z;
    g_8_5_z = g_8_1_z+1;
    g_8_6_z = g_8_5_z;
    g_8_7_z = g_8_5_z;
    g_8_8_z = g_8_5_z;
end
%1 4 5
Jac_6 = reshape(vent_jac(:,1),im_sz);
%1 3 5
Jac_7 = reshape(vent_jac(:,2),im_sz);
%2 4 5
Jac_2 = reshape(vent_jac(:,3),im_sz);
%2 3 5
Jac_3 = reshape(vent_jac(:,4),im_sz);
%1 4 6
Jac_5 = reshape(vent_jac(:,5),im_sz);
%2 3 6
Jac_4 = reshape(vent_jac(:,6),im_sz);
%2 4 6
Jac_1 = reshape(vent_jac(:,7),im_sz);
%1 3 6
Jac_8 = reshape(vent_jac(:,8),im_sz);

coef = 0.5*0.5*0.5;
jac_g_1_1 = Jac_1(g_1_1_x,g_1_1_y,g_1_1_z); 
jac_g_1_2 = Jac_2(g_1_2_x,g_1_2_y,g_1_2_z);
jac_g_1_3 = Jac_3(g_1_3_x,g_1_3_y,g_1_3_z);
jac_g_1_4 = Jac_4(g_1_4_x,g_1_4_y,g_1_4_z);
jac_g_1_5 = Jac_5(g_1_5_x,g_1_5_y,g_1_5_z);
jac_g_1_6 = Jac_6(g_1_6_x,g_1_6_y,g_1_6_z);
jac_g_1_7 = Jac_7(g_1_7_x,g_1_7_y,g_1_7_z);
jac_g_1_8 = Jac_8(g_1_8_x,g_1_8_y,g_1_8_z);
Jac_group_1 = (jac_g_1_1+jac_g_1_2+jac_g_1_3+jac_g_1_4+jac_g_1_5+jac_g_1_6+jac_g_1_7+jac_g_1_8)/8;

jac_g_2_1 = Jac_1(g_2_1_x,g_2_1_y,g_2_1_z); 
jac_g_2_2 = Jac_2(g_2_2_x,g_2_2_y,g_2_2_z);
jac_g_2_3 = Jac_3(g_2_3_x,g_2_3_y,g_2_3_z);
jac_g_2_4 = Jac_4(g_2_4_x,g_2_4_y,g_2_4_z);
jac_g_2_5 = Jac_5(g_2_5_x,g_2_5_y,g_2_5_z);
jac_g_2_6 = Jac_6(g_2_6_x,g_2_6_y,g_2_6_z);
jac_g_2_7 = Jac_7(g_2_7_x,g_2_7_y,g_2_7_z);
jac_g_2_8 = Jac_8(g_2_8_x,g_2_8_y,g_2_8_z);
Jac_group_2 = (jac_g_2_1+jac_g_2_2+jac_g_2_3+jac_g_2_4+jac_g_2_5+jac_g_2_6+jac_g_2_7+jac_g_2_8)/8;

jac_g_3_1 = Jac_1(g_3_1_x,g_3_1_y,g_3_1_z); 
jac_g_3_2 = Jac_2(g_3_2_x,g_3_2_y,g_3_2_z);
jac_g_3_3 = Jac_3(g_3_3_x,g_3_3_y,g_3_3_z);
jac_g_3_4 = Jac_4(g_3_4_x,g_3_4_y,g_3_4_z);
jac_g_3_5 = Jac_5(g_3_5_x,g_3_5_y,g_3_5_z);
jac_g_3_6 = Jac_6(g_3_6_x,g_3_6_y,g_3_6_z);
jac_g_3_7 = Jac_7(g_3_7_x,g_3_7_y,g_3_7_z);
jac_g_3_8 = Jac_8(g_3_8_x,g_3_8_y,g_3_8_z);
Jac_group_3 = (jac_g_3_1+jac_g_3_2+jac_g_3_3+jac_g_3_4+jac_g_3_5+jac_g_3_6+jac_g_3_7+jac_g_3_8)/8;


jac_g_4_1 = Jac_1(g_4_1_x,g_4_1_y,g_4_1_z); 
jac_g_4_2 = Jac_2(g_4_2_x,g_4_2_y,g_4_2_z);
jac_g_4_3 = Jac_3(g_4_3_x,g_4_3_y,g_4_3_z);
jac_g_4_4 = Jac_4(g_4_4_x,g_4_4_y,g_4_4_z);
jac_g_4_5 = Jac_5(g_4_5_x,g_4_5_y,g_4_5_z);
jac_g_4_6 = Jac_6(g_4_6_x,g_4_6_y,g_4_6_z);
jac_g_4_7 = Jac_7(g_4_7_x,g_4_7_y,g_4_7_z);
jac_g_4_8 = Jac_8(g_4_8_x,g_4_8_y,g_4_8_z);
Jac_group_4 = (jac_g_4_1+jac_g_4_2+jac_g_4_3+jac_g_4_4+jac_g_4_5+jac_g_4_6+jac_g_4_7+jac_g_4_8)/8;


jac_g_5_1 = Jac_1(g_5_1_x,g_5_1_y,g_5_1_z); 
jac_g_5_2 = Jac_2(g_5_2_x,g_5_2_y,g_5_2_z);
jac_g_5_3 = Jac_3(g_5_3_x,g_5_3_y,g_5_3_z);
jac_g_5_4 = Jac_4(g_5_4_x,g_5_4_y,g_5_4_z);
jac_g_5_5 = Jac_5(g_5_5_x,g_5_5_y,g_5_5_z);
jac_g_5_6 = Jac_6(g_5_6_x,g_5_6_y,g_5_6_z);
jac_g_5_7 = Jac_7(g_5_7_x,g_5_7_y,g_5_7_z);
jac_g_5_8 = Jac_8(g_5_8_x,g_5_8_y,g_5_8_z);
Jac_group_5 = (jac_g_5_1+jac_g_5_2+jac_g_5_3+jac_g_5_4+jac_g_5_5+jac_g_5_6+jac_g_5_7+jac_g_5_8)/8;

jac_g_6_1 = Jac_1(g_6_1_x,g_6_1_y,g_6_1_z); 
jac_g_6_2 = Jac_2(g_6_2_x,g_6_2_y,g_6_2_z);
jac_g_6_3 = Jac_3(g_6_3_x,g_6_3_y,g_6_3_z);
jac_g_6_4 = Jac_4(g_6_4_x,g_6_4_y,g_6_4_z);
jac_g_6_5 = Jac_5(g_6_5_x,g_6_5_y,g_6_5_z);
jac_g_6_6 = Jac_6(g_6_6_x,g_6_6_y,g_6_6_z);
jac_g_6_7 = Jac_7(g_6_7_x,g_6_7_y,g_6_7_z);
jac_g_6_8 = Jac_8(g_6_8_x,g_6_8_y,g_6_8_z);
Jac_group_6 = (jac_g_6_1+jac_g_6_2+jac_g_6_3+jac_g_6_4+jac_g_6_5+jac_g_6_6+jac_g_6_7+jac_g_6_8)/8;

jac_g_7_1 = Jac_1(g_7_1_x,g_7_1_y,g_7_1_z); 
jac_g_7_2 = Jac_2(g_7_2_x,g_7_2_y,g_7_2_z);
jac_g_7_3 = Jac_3(g_7_3_x,g_7_3_y,g_7_3_z);
jac_g_7_4 = Jac_4(g_7_4_x,g_7_4_y,g_7_4_z);
jac_g_7_5 = Jac_5(g_7_5_x,g_7_5_y,g_7_5_z);
jac_g_7_6 = Jac_6(g_7_6_x,g_7_6_y,g_7_6_z);
jac_g_7_7 = Jac_7(g_7_7_x,g_7_7_y,g_7_7_z);
jac_g_7_8 = Jac_8(g_7_8_x,g_7_8_y,g_7_8_z);
Jac_group_7 = (jac_g_7_1+jac_g_7_2+jac_g_7_3+jac_g_7_4+jac_g_7_5+jac_g_7_6+jac_g_7_7+jac_g_7_8)/8;

jac_g_8_1 = Jac_1(g_8_1_x,g_8_1_y,g_8_1_z); 
jac_g_8_2 = Jac_2(g_8_2_x,g_8_2_y,g_8_2_z);
jac_g_8_3 = Jac_3(g_8_3_x,g_8_3_y,g_8_3_z);
jac_g_8_4 = Jac_4(g_8_4_x,g_8_4_y,g_8_4_z);
jac_g_8_5 = Jac_5(g_8_5_x,g_8_5_y,g_8_5_z);
jac_g_8_6 = Jac_6(g_8_6_x,g_8_6_y,g_8_6_z);
jac_g_8_7 = Jac_7(g_8_7_x,g_8_7_y,g_8_7_z);
jac_g_8_8 = Jac_8(g_8_8_x,g_8_8_y,g_8_8_z);
Jac_group_8 = (jac_g_8_1+jac_g_8_2+jac_g_8_3+jac_g_8_4+jac_g_8_5+jac_g_8_6+jac_g_8_7+jac_g_8_8)/8;

vent_img = zeros(im_sz(1)-1,img_sz(2)-1,im_sz(3)-1);
vent_img(g_1_1_x,g_1_1_y,g_1_1_z) = Jac_group_1;
vent_img(g_2_1_x,g_2_1_y,g_2_1_z) = Jac_group_2;
vent_img(g_3_1_x,g_3_1_y,g_3_1_z) = Jac_group_3;
vent_img(g_4_1_x,g_4_1_y,g_4_1_z) = Jac_group_4;
vent_img(g_5_1_x,g_5_1_y,g_5_1_z) = Jac_group_5;
vent_img(g_6_1_x,g_6_1_y,g_6_1_z) = Jac_group_6;
vent_img(g_7_1_x,g_7_1_y,g_7_1_z) = Jac_group_7;
vent_img(g_8_1_x,g_8_1_y,g_8_1_z) = Jac_group_8;




