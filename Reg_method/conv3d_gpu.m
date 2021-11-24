function [gy_filter,tar_sum] = conv3d_gpu(img,filter,m_filter,groups)

use_cmd = 0;
if use_cmd ~= 1
    
    [N,channel,r,c,s] = size(img);
    [N,channel,fr,fc,fs] = size(filter);
    img_1 = [];
    filter_1 = [];
    for i = 1:channel
        img_cur = reshape(img(:,i,:,:,:),1,N*r*c*s);
        img_1 = cat(2,img_1,img_cur);
        img_1 = single(img_1);
        filter_cur = reshape(filter(:,i,:,:,:),1,N*fr*fc*fs);
        filter_1 = cat(2,filter_1,filter_cur);
        filter_1 = single(filter_1);
    end
    output = py.conv_xp.conv3D_accelerate(img_1,filter_1,groups,channel,r,c,s,fr,fc,fs);
%     output_1 = py.conv_1.conv3D(img_1,filter_1,groups,channel,r,c,s,fr,fc,fs);
    
    gy_filter = cat(1, double(output{1}), double(output{2}), double(output{3}),...
                double(output{4}), double(output{5}), double(output{6}));
    tar_sum = cat(1, double(output{7}), double(output{8}), double(output{9}),...
                double(output{10}), double(output{11}), double(output{12}));
    gy_filter = squeeze(mean(gy_filter,1));
    tar_sum = squeeze(mean(tar_sum,1));
else
    save inter.mat img filter m_filter groups -v7
    path = cd;
    opean_path = ['CD',' ',path];
    activate_env = ['conda activate Conv3D_GPU'];
    file = ['python conv3d_gpu_cuda.py'];
    tic
    commd = [opean_path '&&' activate_env '&&' file];
    state = system(commd);
    toc
    if state ~= 0
        error(['cmd error.']); 
    end
    load('results.mat');
    conv = squeeze(array);
    gy_filter = conv(7:12,:,:,:,:);
    tar_sum = conv(1:6,:,:,:,:);
    gy_filter = squeeze(mean(gy_filter,1));
    tar_sum = squeeze(mean(tar_sum,1));
end
end

