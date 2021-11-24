Target = randn(100,6,11,11,11);
filter = randn(100,6,7,7,7);
Target_r = reshape(Target,1,600,11,11,11);
group = 100;
path = cd;
opean_path = ['CD',' ',path];
activate_env = ['conda activate Conv3D_GPU'];
file = ['python conv3d_gpu_TEST1.py'];
system(opean);
commd = [opean_path '&&' activate_env '&&' file];
state = system(commd);

save inter.mat Target_r group filter -v7
