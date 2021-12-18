clear;
%This code is based on pre_process_train_v5.m , generate
%UII_preprocess include nor resize cut .etc
addpath(genpath('../mutils/My/'));
addpath(genpath('../ptv'));
addpath(genpath('scr/'))
data_path = 'D:/CTVI_data/Step1_pre_processing/';
% data_path = 'D:/LiverDataset/UII_liver_all/';
%%
% preprocess
count_all_sz = zeros(46,3);
spc_count = zeros(46,3);
for n_index = 1:46

    info = niftiinfo([data_path,'Sub_',num2str(n_index),'_CT_avg_mask.nii']);

    count_all_sz(n_index,:) =info.ImageSize;
       
end
min_sz = min(count_all_sz);
max_sz = max(count_all_sz);
m_sz = mean(count_all_sz);
final_sz = ceil(ceil(m_sz/2)/16)*16;
