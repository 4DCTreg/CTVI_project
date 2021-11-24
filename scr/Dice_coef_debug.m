load('debug_dice.mat');
Compare_vent_jac = reshape(vent_jac,[GT_img_sz(1)*GT_img_sz(2)*GT_img_sz(3),1]);
Compare_vent_GT = reshape(vent_GT,[GT_img_sz(1)*GT_img_sz(2)*GT_img_sz(3),1]);
vent_index = Compare_vent_jac;
gt_index = Compare_vent_GT;
vent_mask = find(Compare_vent_jac~=0);
gt_mask = find(Compare_vent_GT~=0);
Compare_vent_jac(find(Compare_vent_jac==0)) = [];
Compare_vent_jac = img_thr_v1(Compare_vent_jac, min(Compare_vent_jac), max(Compare_vent_jac), 15000);
vent_index(vent_mask) = Compare_vent_jac;
vent_index_sort = sort(Compare_vent_jac);
Compare_vent_GT(find(Compare_vent_GT==0)) = [];
Compare_vent_GT = img_thr_v1(Compare_vent_GT, min(Compare_vent_GT), max(Compare_vent_GT) , 15000);
gt_index(gt_mask) = Compare_vent_GT;
cp_sz = size(Compare_vent_GT,1);
gt_index_sort = sort(Compare_vent_GT);
vent_0 = vent_index_sort(1);
vent_2 = vent_index_sort(round(0.2*cp_sz));
vent_4 = vent_index_sort(round(0.4*cp_sz));
vent_6 = vent_index_sort(round(0.6*cp_sz));
vent_8 = vent_index_sort(round(0.8*cp_sz));
vent_10 = vent_index_sort(end);

gt_0 = gt_index_sort(1);
gt_2 = gt_index_sort(round(0.2*cp_sz));
gt_4 = gt_index_sort(round(0.4*cp_sz));
gt_6 = gt_index_sort(round(0.6*cp_sz));
gt_8 = gt_index_sort(round(0.8*cp_sz));
gt_10 = gt_index_sort(end);
v_D1 = find(vent_index>=vent_0&vent_index<vent_2);
v_D2 = find(vent_index>=vent_2&vent_index<vent_4);
v_D3 = find(vent_index>=vent_4&vent_index<vent_6);
v_D4 = find(vent_index>=vent_6&vent_index<vent_8);
v_D5 = find(vent_index>=vent_8&vent_index<vent_10);

g_D1 = find(gt_index>=gt_0&gt_index<gt_2);
g_D2 = find(gt_index>=gt_2&gt_index<gt_4);
g_D3 = find(gt_index>=gt_4&gt_index<gt_6);
g_D4 = find(gt_index>=gt_6&gt_index<gt_8);
g_D5 = find(gt_index>=gt_8&gt_index<gt_10);

intersect1 = intersect(v_D1,g_D1);
intersect2 = intersect(v_D2,g_D2);
intersect3 = intersect(v_D3,g_D3);
intersect4 = intersect(v_D4,g_D4);
intersect5 = intersect(v_D5,g_D5);
D1 = 2*size(intersect1,1)/(size(v_D1,1)+size(g_D1,1));
D2 = 2*size(intersect2,1)/(size(v_D2,1)+size(g_D2,1));
D3 = 2*size(intersect3,1)/(size(v_D3,1)+size(g_D3,1));
D4 = 2*size(intersect4,1)/(size(v_D4,1)+size(g_D4,1));
D5 = 2*size(intersect5,1)/(size(v_D5,1)+size(g_D5,1));
