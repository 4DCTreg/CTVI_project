clc;
clear;

fid = fopen('VAMPIRE_test_all.txt');
tline = fgetl(fid);
list_cell={};
list_cell = [list_cell;tline];
while ischar(tline)
    tline = fgetl(fid);
    list_cell = [list_cell;tline];
end
fclose(fid);


results_path = 'D:/CTVI_data/Gan_results';
% info_Ap = niftiinfo([data_path,'image_Ap',num2str(n_index),'.nii']);
for n_index =  1 : 7
    GT_name_cur = list_cell(n_index*6); 
    GT_img_name = GT_name_cur{1}(34:end);
    results_name_cur = [results_path, GT_img_name];
    info_cur = niftiinfo(results_name_cur);
    img = niftiread(results_name_cur);
    info_cur.raw.qoffset_x = 0;
    info_cur.raw.qoffset_y = 0;
    info_cur.raw.qoffset_z = 0;
    info_cur.Transform.T(4,:) = [1,1,1,1];
    niftiwrite(img,[results_name_cur],info_cur);
    
end