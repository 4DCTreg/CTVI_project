load('results_v1.mat');
count = 0;
for i = 1:25
    cur = landmarks_SIFT{1};
    cur_size = size(cur,1);
    count = count+cur_size;
end