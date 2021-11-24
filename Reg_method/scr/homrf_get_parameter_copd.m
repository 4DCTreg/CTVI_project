function parameter=homrf_get_parameter_copd(idx,useTop,spc)
switch idx
    case 1
        parameter.nlevel=7;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1; 1 1 1; 1 1 1; 2 2 1; 2 2 1; 2 3 2;2 3 2];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 1;2 2 2;4 4 2;6 6 3; 8 8 4;10 10 5];
        parameter.x_tot=201;
        parameter.y_tot=201;
        parameter.z_tot=201;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 9;11 11 11];
        parameter.presmooth=0;
        parameter.spc=spc;
        parameter.smooth_co=[1,100,20,15,3,3,2];
        parameter.Tcoe_n=[20,5,5,5,5,5,1];
        parameter.end=2;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)*(12/5) ,spc(1)*spc(2)*spc(3)*(12/5),  spc(1)*spc(2)*spc(3)*(12/6),  spc(1)*spc(2)*spc(3)*(12/8),  spc(1)*spc(2)*spc(3)*(12/10),  spc(1)*spc(2)*spc(3)/(12/12)^3;];
        parameter.resize=0;
        parameter.dist_co=100;

    case 2
        parameter.nlevel=7;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;2 2 2];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 1;4 4 2;6 6 3;8 8 4; 10 10 5;12 12 6];
        parameter.x_tot=171;
        parameter.y_tot=171;
        parameter.z_tot=171;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 7 5;7 7 7;7 11 9;7 13 9;9 15 13;7 13 11];
        parameter.presmooth=0;
        parameter.spc=spc;
        parameter.smooth_co=[1,50,20,5,3,2,1];
        parameter.Tcoe_n=[20,0.005,0.001,0.001,0.05,0.05,0.1];
        parameter.end=2;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)*(12/5) ,spc(1)*spc(2)*spc(3)*(12/5),  spc(1)*spc(2)*spc(3)*(12/6),  spc(1)*spc(2)*spc(3)*(12/8),  spc(1)*spc(2)*spc(3)*(12/10),  spc(1)*spc(2)*spc(3)/(12/12)^3;];
        parameter.resize=0;
        parameter.dist_co=100;
        
    case 3
        parameter.nlevel=7;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;1 1 1];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 1;3 3 2;4 4 3;5 5 4; 6 6 5;7 7 6];
        parameter.x_tot=201;
        parameter.y_tot=201;
        parameter.z_tot=201;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 9;5 5 3];
        parameter.presmooth=0;
        parameter.spc=spc;
        parameter.smooth_co=[1,100,20,15,3,3,2];
        parameter.Tcoe_n=[20,0.005,0.001,0.001,0.05,0.05,0.1];
        parameter.end=2;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)*(12/5) ,spc(1)*spc(2)*spc(3)*(12/5),  spc(1)*spc(2)*spc(3)*(12/6),  spc(1)*spc(2)*spc(3)*(12/8),  spc(1)*spc(2)*spc(3)*(12/10),  spc(1)*spc(2)*spc(3)/(12/12)^3;];
        parameter.resize=0;
        parameter.dist_co=100;
        
    case 4
        parameter.nlevel=7;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;1 1 1];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 1;3 3 2;4 4 3;5 5 4; 6 6 5;7 7 6];
        parameter.x_tot=201;
        parameter.y_tot=201;
        parameter.z_tot=201;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 9;5 5 3];
        parameter.presmooth=0;
        parameter.spc=spc;
        parameter.smooth_co=[1,100,20,15,3,3,2];
        parameter.Tcoe_n=[20,0.005,0.001,0.001,0.05,0.05,0.1];
        parameter.end=2;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)*(12/5) ,spc(1)*spc(2)*spc(3)*(12/5),  spc(1)*spc(2)*spc(3)*(12/6),  spc(1)*spc(2)*spc(3)*(12/8),  spc(1)*spc(2)*spc(3)*(12/10),  spc(1)*spc(2)*spc(3)/(12/12)^3;];
        parameter.resize=0;
        parameter.dist_co=100;
        
        
    case 5
        parameter.nlevel=7;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;1 1 1];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 1;3 3 2;4 4 3;5 5 4; 6 6 5;7 7 6];
        parameter.x_tot=201;
        parameter.y_tot=201;
        parameter.z_tot=201;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 9;5 5 3];
        parameter.presmooth=0;
        parameter.spc=spc;
        parameter.smooth_co=[1,100,20,15,3,3,2];
        parameter.Tcoe_n=[20,0.005,0.001,0.001,0.05,0.05,0.1];
        parameter.end=2;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)*(12/5) ,spc(1)*spc(2)*spc(3)*(12/5),  spc(1)*spc(2)*spc(3)*(12/6),  spc(1)*spc(2)*spc(3)*(12/8),  spc(1)*spc(2)*spc(3)*(12/10),  spc(1)*spc(2)*spc(3)/(12/12)^3;];
        parameter.resize=0;
        parameter.dist_co=100;
       
    case 6
        parameter.nlevel=7;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;1 1 1];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 1;3 3 2;4 4 3;5 5 4; 6 6 5;7 7 6];
        parameter.x_tot=201;
        parameter.y_tot=201;
        parameter.z_tot=201;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 9;5 5 3];
        parameter.presmooth=0;
        parameter.spc=spc;
        parameter.smooth_co=[1,100,20,15,3,3,2];
        parameter.Tcoe_n=[20,0.005,0.001,0.001,0.05,0.05,0.1];
        parameter.end=2;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)*(12/5) ,spc(1)*spc(2)*spc(3)*(12/5),  spc(1)*spc(2)*spc(3)*(12/6),  spc(1)*spc(2)*spc(3)*(12/8),  spc(1)*spc(2)*spc(3)*(12/10),  spc(1)*spc(2)*spc(3)/(12/12)^3;];
        parameter.resize=0;
        parameter.dist_co=100;
       
    case 7
        parameter.nlevel=7;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;1 1 1];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 1;3 3 2;4 4 3;5 5 4; 6 6 5;7 7 6];
        parameter.x_tot=201;
        parameter.y_tot=201;
        parameter.z_tot=201;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 9;5 5 3];
        parameter.presmooth=0;
        parameter.spc=spc;
        parameter.smooth_co=[1,100,20,15,3,3,2];
        parameter.Tcoe_n=[20,0.005,0.001,0.001,0.05,0.05,0.1];
        parameter.end=2;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)*(12/5) ,spc(1)*spc(2)*spc(3)*(12/5),  spc(1)*spc(2)*spc(3)*(12/6),  spc(1)*spc(2)*spc(3)*(12/8),  spc(1)*spc(2)*spc(3)*(12/10),  spc(1)*spc(2)*spc(3)/(12/12)^3;];
        parameter.resize=0;
        parameter.dist_co=100;
    case 8
        parameter.nlevel=7;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;1 1 1];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 1;3 3 2;4 4 3;5 5 4; 6 6 5;7 7 6];
        parameter.x_tot=201;
        parameter.y_tot=201;
        parameter.z_tot=201;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 9;5 5 3];
        parameter.presmooth=0;
        parameter.spc=spc;
        parameter.smooth_co=[1,100,20,15,3,3,2];
        parameter.Tcoe_n=[20,0.005,0.001,0.001,0.05,0.05,0.1];
        parameter.end=2;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)*(12/5) ,spc(1)*spc(2)*spc(3)*(12/5),  spc(1)*spc(2)*spc(3)*(12/6),  spc(1)*spc(2)*spc(3)*(12/8),  spc(1)*spc(2)*spc(3)*(12/10),  spc(1)*spc(2)*spc(3)/(12/12)^3;];
        parameter.resize=0;
        parameter.dist_co=100;
    case 9
        parameter.nlevel=7;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;1 1 1];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 1;3 3 2;4 4 3;5 5 4; 6 6 5;7 7 6];
        parameter.x_tot=201;
        parameter.y_tot=201;
        parameter.z_tot=201;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 9;5 5 3];
        parameter.presmooth=0;
        parameter.spc=spc;
        parameter.smooth_co=[1,100,20,15,3,3,2];
        parameter.Tcoe_n=[20,0.005,0.001,0.001,0.05,0.05,0.1];
        parameter.end=2;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)*(12/5) ,spc(1)*spc(2)*spc(3)*(12/5),  spc(1)*spc(2)*spc(3)*(12/6),  spc(1)*spc(2)*spc(3)*(12/8),  spc(1)*spc(2)*spc(3)*(12/10),  spc(1)*spc(2)*spc(3)/(12/12)^3;];
        parameter.resize=0;
        parameter.dist_co=100;
       
    case 10
        parameter.nlevel=7;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;1 1 1];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 1;3 3 2;4 4 3;5 5 4; 6 6 5;7 7 6];
        parameter.x_tot=201;
        parameter.y_tot=201;
        parameter.z_tot=201;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 9;5 5 3];
        parameter.presmooth=0;
        parameter.spc=spc;
        parameter.smooth_co=[1,100,20,15,3,3,2];
        parameter.Tcoe_n=[20,0.005,0.001,0.001,0.05,0.05,0.1];
        parameter.end=2;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)*(12/5) ,spc(1)*spc(2)*spc(3)*(12/5),  spc(1)*spc(2)*spc(3)*(12/6),  spc(1)*spc(2)*spc(3)*(12/8),  spc(1)*spc(2)*spc(3)*(12/10),  spc(1)*spc(2)*spc(3)/(12/12)^3;];
        parameter.resize=0;
        parameter.dist_co=100;
end
end