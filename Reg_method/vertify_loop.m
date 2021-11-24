clc;
clear;
for i = 1:1000
target = randn(7,7,3,7);
target_add = randn(1,10);
target_add_1 = randn(1,10);
r1 = py.fun_test.add(target_add,target_add_1);
r2 = double(r1{1});
disp(['iter=' num2str(i) ]);
end



