function vent_img = jac_interp(im_sz,vent_jac)
% vent_img(find(vent_img)<0) = 0;
% 3D linear interpolation
% divide into 8 groups, each groups has 8 points
% groups
%   1—4
%  /  /
% 2—3
%   5—8
%  /  /
% 6—7
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

vent_img = zeros(im_sz(1)-1,im_sz(2)-1,im_sz(3)-1);
vent_img(g_1_1_x,g_1_1_y,g_1_1_z) = Jac_group_1;
vent_img(g_2_1_x,g_2_1_y,g_2_1_z) = Jac_group_2;
vent_img(g_3_1_x,g_3_1_y,g_3_1_z) = Jac_group_3;
vent_img(g_4_1_x,g_4_1_y,g_4_1_z) = Jac_group_4;
vent_img(g_5_1_x,g_5_1_y,g_5_1_z) = Jac_group_5;
vent_img(g_6_1_x,g_6_1_y,g_6_1_z) = Jac_group_6;
vent_img(g_7_1_x,g_7_1_y,g_7_1_z) = Jac_group_7;
vent_img(g_8_1_x,g_8_1_y,g_8_1_z) = Jac_group_8;
vent_img = volresize(vent_img,im_sz,1);
end