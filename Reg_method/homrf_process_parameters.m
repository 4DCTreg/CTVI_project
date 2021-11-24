function [nlevel,k_down,x_tot, y_tot, z_tot, useTop, grid_space, metric, labels_search,...
    presmooth, quant, smooth_co_pre, endl, Tcoe_n, top_co, spc, resize, dist_co] = homrf_process_parameters(parameters)
nlevel = getoptions(parameters, 'nlevel',3);
k_down = getoptions(parameters, 'k_down',0.8);
x_tot = getoptions(parameters, 'x_tot',103);
y_tot = getoptions(parameters, 'y_tot',103);
z_tot= getoptions(parameters, 'z_tot',103);
useTop= getoptions(parameters, 'useTop',0);
grid_space= getoptions(parameters, 'grid_space',[ 6 6 3; 9 9 4; 12 12 5; ]);
metric= getoptions(parameters, 'metric','MIND');
labels_search= getoptions(parameters, 'labels',[ 7 7 7;9 9 9; 11 11 11;]);
presmooth= getoptions(parameters, 'presmooth', []);
quant= getoptions(parameters, 'quant',[ 1 1 2; 1 1 2; 2 2 3]);
smooth_co_pre= getoptions(parameters, 'smooth_co',0.1);
endl = getoptions(parameters, 'endl',1);
Tcoe_n= getoptions(parameters, 'Tcoe_n',0.8);
top_co= getoptions(parameters, 'top_co',0.8);
spc= getoptions(parameters, 'spc',[]);
resize= getoptions(parameters, 'resize',1);
dist_co= getoptions(parameters, 'dist_co',1);
end