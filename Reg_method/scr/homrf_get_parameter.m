function parameter=homrf_get_parameter(idx,useTop,spc)
switch idx
    case 1
        parameter.nlevel=1;
        parameter.k_down=0.8;
        parameter.quant=[ 1 1 1];
        parameter.useTop=useTop;
        parameter.grid_space=[ 8 8 8];
        parameter.x_tot=33;
        parameter.y_tot=33;
        parameter.z_tot=33;
        parameter.metric='MIND';
        parameter.labels=[11 11 11;];
        parameter.presmooth=0;
        parameter.spc=spc;
        if useTop==1
            parameter.smooth_co=0.05;
        else
            parameter.smooth_co=0.05;
        end
        parameter.end=1;
        parameter.top_co= 1;
        parameter.Tcoe_n=5;
        parameter.resize=1;
        parameter.dist_co=1;
    case 2
        parameter.nlevel=6;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1;1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 2;4 4 4;4 4 4; 6 6 6; 8 8 8];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 13;];
        parameter.presmooth=0;
        parameter.spc=spc;
        if useTop==1
            parameter.smooth_co=0.05;
        else
            parameter.smooth_co=0.05;
        end
        parameter.end=6;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3   ,spc(1)*spc(2)*spc(3)/(10/2)^3,  spc(1)*spc(2)*spc(3)/(10/2)^3,    spc(1)*spc(2)*spc(3)/(10/6)^3,  spc(1)*spc(2)*spc(3)/(10/8)^3,    spc(1)*spc(2)*spc(3)/1^3;];
        parameter.top_co=[1   ,0.008,  0.06,  0.125,  0.42, 1;];
        parameter.Tcoe_n=5;
        parameter.resize=1;
        parameter.dist_co=1;
%         parameter.nlevel=6;
%         parameter.k_down=0.8;
%         parameter.quant=[1 1 1;1 1 1; 1 1 1; 1 1 1; 1 1 2; 1 1 1];
%         parameter.useTop=useTop;
%         parameter.grid_space=[2 2 2;2 2 2;4 4 4;6 6 6;8 8 8; 8 8 8];
%         parameter.x_tot=33;
%         parameter.y_tot=33;
%         parameter.z_tot=33;
%         parameter.metric='MIND';
%         parameter.labels=[3 3 3;3 3 3;3 3 3;7 7 7;9 9 9;9 9 9;];
%         parameter.presmooth=1;
%         parameter.spc=spc;
%         parameter.smooth_co=[1,20,5,10/6,10/8,10/10];
%         parameter.Tcoe_n=[20,0.05,0.05,0.05,0.05,0.05];
%         parameter.end=2;
%         parameter.resize=0;
%         parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)/100,  spc(1)*spc(2)*spc(3)/64,  spc(1)*spc(2)*spc(3)/(10/6)^3,  spc(1)*spc(2)*spc(3)/(10/8)^3,  spc(1)*spc(2)*spc(3)/1^3;];
%         parameter.dist_co=30.76;
    case 3
%         parameter.nlevel=6;
%         parameter.k_down=0.8;
%         parameter.quant=[1 1 1;1 1 1; 1 1 1; 1 1 2; 1 1 3; 2 2 4];
%         parameter.useTop=useTop;
%         parameter.grid_space=[2 2 1;2 2 2;4 4 4 ;6 6 6; 8 8 8; 10 10 10];
%         parameter.x_tot=103;
%         parameter.y_tot=103;
%         parameter.z_tot=103;
%         parameter.metric='MIND';
%         parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 9;9 9 9;];
%         parameter.presmooth=0;
%         parameter.spc=spc;
%         parameter.smooth_co=[1,10,2.8,10/6,10/8,10/10];
%         parameter.end=3;
%         parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3  ,spc(1)*spc(2)*spc(3)/125,  spc(1)*spc(2)*spc(3)/(10/2)^3,   spc(1)*spc(2)*spc(3)/(10/6)^3,  spc(1)*spc(2)*spc(3)/(10/8)^3,   spc(1)*spc(2)*spc(3)/1^3;];
%         parameter.Tcoe_n=[20,5,5,5,5,5];
%         parameter.resize=1;
%         parameter.dist_co=30.76;
        
        
        parameter.nlevel=6;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1;1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 2;4 4 4;4 4 4; 6 6 6; 8 8 8];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 13;];
        parameter.presmooth=0;
        parameter.spc=spc;
        if useTop==1
            parameter.smooth_co=0.05;
        else
            parameter.smooth_co=0.05;
        end
        parameter.end=6;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3   ,spc(1)*spc(2)*spc(3)/(10/2)^3,  spc(1)*spc(2)*spc(3)/(10/2)^3,    spc(1)*spc(2)*spc(3)/(10/6)^3,  spc(1)*spc(2)*spc(3)/(10/8)^3,    spc(1)*spc(2)*spc(3)/1^3;];
        parameter.top_co=[1   ,0.008,  0.06,  0.125,  0.42, 1;];
        parameter.Tcoe_n=5;
        parameter.resize=1;
        parameter.dist_co=1;
    case 4
%         parameter.nlevel=6;
%         parameter.k_down=0.8;
%         parameter.quant=[1 1 1;1 1 2; 1 1 2; 1 1 2; 2 2 3; 2 2 4];
%         parameter.useTop=useTop;
%         parameter.grid_space=[2 2 2;2 2 2;4 4 4;6 6 6; 8 8 8; 10 10 10];
%         parameter.x_tot=103;
%         parameter.y_tot=103;
%         parameter.z_tot=103;
%         parameter.metric='MIND';
%         parameter.labels=[3 3 3; 7 7 7; 9 9 9; 11 11 11; 11 9 9; 11 11 11;];
%         parameter.presmooth=0;
%         parameter.spc=spc;
%         parameter.smooth_co=[1,18,15,10/6,10/8,10/10];
%         parameter.Tcoe_n=[20,5,5,5,5,5];
%         parameter.end=2;
%         parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)/(10/2)^3,  spc(1)*spc(2)*spc(3)/64,  spc(1)*spc(2)*spc(3)/(10/6)^3,  spc(1)*spc(2)*spc(3)/(10/8)^3,  spc(1)*spc(2)*spc(3)/1^3;];
%         parameter.resize=1;
%         parameter.dist_co=30.76;
        
        parameter.nlevel=6;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1;1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 2;4 4 4;4 4 4; 6 6 6; 8 8 8];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 13;];
        parameter.presmooth=0;
        parameter.spc=spc;
        if useTop==1
            parameter.smooth_co=0.05;
        else
            parameter.smooth_co=0.05;
        end
        parameter.end=6;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3   ,spc(1)*spc(2)*spc(3)/(10/2)^3,  spc(1)*spc(2)*spc(3)/(10/2)^3,    spc(1)*spc(2)*spc(3)/(10/6)^3,  spc(1)*spc(2)*spc(3)/(10/8)^3,    spc(1)*spc(2)*spc(3)/1^3;];
        parameter.top_co=[1   ,0.008,  0.06,  0.125,  0.42, 1;];
        parameter.Tcoe_n=5;
        parameter.resize=1;
        parameter.dist_co=1;
    case 5
                parameter.nlevel=6;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1;1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 2;4 4 4;4 4 4; 6 6 6; 8 8 8];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;13 13 13;];
        parameter.presmooth=0;
        parameter.spc=spc;
        if useTop==1
            parameter.smooth_co=0.05;
        else
            parameter.smooth_co=0.05;
        end
        parameter.end=6;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3   ,spc(1)*spc(2)*spc(3)/(10/2)^3,  spc(1)*spc(2)*spc(3)/(10/2)^3,    spc(1)*spc(2)*spc(3)/(10/6)^3,  spc(1)*spc(2)*spc(3)/(10/8)^3,    spc(1)*spc(2)*spc(3)/1^3;];
        parameter.top_co=[1   ,0.008,  0.06,  0.125,  0.42, 1;];
        parameter.Tcoe_n=5;
        parameter.resize=1;
        parameter.dist_co=1;
%         parameter.nlevel=4;
%         parameter.k_down=0.8;
%         parameter.quant=[ 1 1 1;  1 1 1; 1 1 1;  1 1 1];
%         parameter.useTop=useTop;
%         parameter.grid_space=[ 4 4 4; 6 6 6; 8 8 8; 8 8 8];
%         parameter.x_tot=103;
%         parameter.y_tot=103;
%         parameter.z_tot=103;
%         parameter.metric='MIND';
%         parameter.labels=[5 5 5;7 7 7;9 9 9;11 11 11];
%         parameter.presmooth=0;
%         parameter.spc=spc;
% %         parameter.smooth_co=[1, 1, 2.5, 0.31, 0.09, 0.04, 0.02];
%         parameter.smooth_co=0.02;
% %         parameter.Tcoe_n=[20,0.005,0.001,0.001,0.05,0.05,0.1];
%         parameter.top_co=[5,  5,...
%            5,  5];
%         parameter.Tcoe_n=1;
%         parameter.end=1;
% %         parameter.top_co=1;
%         parameter.resize=1;
%         parameter.dist_co=1;
    case 6
%         parameter.nlevel=7;
%         parameter.k_down=0.8;
%         parameter.quant=[1 1 2; 1 1 2; 1 1 2; 1 1 2; 1 1 2; 1 1 2;2 2 4];
%         parameter.useTop=useTop;
%         parameter.grid_space=[1 1 1;2 2 2;4 4 4;6 6 6;8 8 8; 10 10 10;12 12 12];
%         parameter.x_tot=101;
%         parameter.y_tot=101;
%         parameter.z_tot=101;
%         parameter.metric='MIND';
%         parameter.labels=[5 5 5;5 5 5;7 7 7;7 7 7;11 11 11;13 13 13;13 13 13];
%         parameter.presmooth=0;
%         parameter.spc=spc;
%         parameter.smooth_co=[100,60,10,5,5,3,1];
%         parameter.Tcoe_n=[20,0.005,0.001,0.001,0.05,0.05,0.01];
%         parameter.end=1;
%         parameter.top_co=[1,1,1,1,1,1,1];
%         parameter.resize=1;
%         parameter.dist_co=100;
        
        parameter.nlevel=6;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1;1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 2;4 4 4;4 4 4; 6 6 6; 8 8 8];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 13;];
        parameter.presmooth=0;
        parameter.spc=spc;
        if useTop==1
            parameter.smooth_co=0.05;
        else
            parameter.smooth_co=0.05;
        end
        parameter.end=6;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3   ,spc(1)*spc(2)*spc(3)/(10/2)^3,  spc(1)*spc(2)*spc(3)/(10/2)^3,    spc(1)*spc(2)*spc(3)/(10/6)^3,  spc(1)*spc(2)*spc(3)/(10/8)^3,    spc(1)*spc(2)*spc(3)/1^3;];
        parameter.top_co=[1   ,0.008,  0.06,  0.125,  0.42, 1;];
        parameter.Tcoe_n=5;
        parameter.resize=1;
        parameter.dist_co=1;
    case 7
        
        parameter.nlevel=6;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1;1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 2;4 4 4;4 4 4; 6 6 6; 8 8 8];
        parameter.x_tot=203;
        parameter.y_tot=203;
        parameter.z_tot=203;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 13;];
        parameter.presmooth=0;
        parameter.spc=spc;
        if useTop==1
            parameter.smooth_co=0.05;
        else
            parameter.smooth_co=0.05;
        end
        parameter.end=6;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3   ,spc(1)*spc(2)*spc(3)/(10/2)^3,  spc(1)*spc(2)*spc(3)/(10/2)^3,    spc(1)*spc(2)*spc(3)/(10/6)^3,  spc(1)*spc(2)*spc(3)/(10/8)^3,    spc(1)*spc(2)*spc(3)/1^3;];
        parameter.top_co=[1   ,0.008,  0.06,  0.125,  0.42, 1;];
        parameter.Tcoe_n=5;
        parameter.resize=1;
        parameter.dist_co=1;
        
%         parameter.nlevel=7;
%         parameter.k_down=0.8;
%         parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;1 1 2];
%         parameter.useTop=useTop;
%         parameter.grid_space=[1 1 1;2 2 1;4 4 2;6 6 3;8 8 4; 10 10 5;12 12  6];
%         parameter.x_tot=71;
%         parameter.y_tot=71;
%         parameter.z_tot=71;
%         parameter.metric='MIND';
%         parameter.labels=[3 3 3;5 5 5;7 7 7;9 9 9;9 9 9;9 13 13;13 13 13];
%         parameter.presmooth=0;
%         parameter.spc=spc;
%         parameter.smooth_co=[1,45,13,4,2.8,2.8,0.8];
%         parameter.Tcoe_n=[20,0.01,0.01,0.01,0.01,0.01,0.01];
%         parameter.end=2;
%         parameter.top_co=[1,2.4/((12/2)^2*6),2.4/(((12/4)^2)*(6/2)),2.4/(((12/6)^2)*(6/3)),2.4/(((12/8)^2)*(6/4)),2.4/(((12/10)^2)*(6/5)),2.4];
%         parameter.resize=0;
%         parameter.dist_co=100;
    case 8
%         parameter.nlevel=7;
%         parameter.k_down=0.8;
%         parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;2 2 2];
%         parameter.useTop=useTop;
%         parameter.grid_space=[2 2 2;2 2 1;4 4 2;6 6 3;8 8 4; 10 10 5;8 8 8];
%         parameter.x_tot=71;
%         parameter.y_tot=71;
%         parameter.z_tot=71;
%         parameter.metric='MIND';
%         parameter.labels=[3 3 3;5 5 5;7 7 7;7 7 5;7 7 7;7 7 7;13 13 13];
%         parameter.presmooth=0;
%         parameter.spc=spc;
%         parameter.smooth_co=[1,80,10,10,1,1,1];
%         parameter.Tcoe_n=[20,0.001,0.001,0.001,0.05,0.05,0.1];
%         parameter.end=2;
%         parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3  ,spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)/(10/2)^3,   spc(1)*spc(2)*spc(3)/64,  spc(1)*spc(2)*spc(3)/(10/6)^3,   spc(1)*spc(2)*spc(3)/(10/8)^3,  spc(1)*spc(2)*spc(3)/1^3;];
%         parameter.resize=0;
%         parameter.dist_co=100;
        
        
        parameter.nlevel=6;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1;1 1 1; 1 1 1; 1 1 1; 1 1 1; 2 2 2];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 2;4 4 4;4 4 4; 6 6 6; 8 8 8];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;13 13 13;];
        parameter.presmooth=0;
        parameter.spc=spc;
        if useTop==1
            parameter.smooth_co=0.05;
        else
            parameter.smooth_co=0.05;
        end
        parameter.end=6;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3   ,spc(1)*spc(2)*spc(3)/(10/2)^3,  spc(1)*spc(2)*spc(3)/(10/2)^3,    spc(1)*spc(2)*spc(3)/(10/6)^3,  spc(1)*spc(2)*spc(3)/(10/8)^3,    spc(1)*spc(2)*spc(3)/1^3;];
        parameter.top_co=[1   ,0.008,  0.06,  0.125,  0.42, 1;];
        parameter.Tcoe_n=5;
        parameter.resize=1;
        parameter.dist_co=1;
    case 9
        
                parameter.nlevel=6;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1;1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 2;4 4 4;4 4 4; 6 6 6; 8 8 8];
        parameter.x_tot=203;
        parameter.y_tot=203;
        parameter.z_tot=203;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 13;];
        parameter.presmooth=0;
        parameter.spc=spc;
        if useTop==1
            parameter.smooth_co=0.05;
        else
            parameter.smooth_co=0.05;
        end
        parameter.end=6;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3   ,spc(1)*spc(2)*spc(3)/(10/2)^3,  spc(1)*spc(2)*spc(3)/(10/2)^3,    spc(1)*spc(2)*spc(3)/(10/6)^3,  spc(1)*spc(2)*spc(3)/(10/8)^3,    spc(1)*spc(2)*spc(3)/1^3;];
        parameter.top_co=[1   ,0.008,  0.06,  0.125,  0.42, 1;];
        parameter.Tcoe_n=5;
        parameter.resize=1;
        parameter.dist_co=1;
        
%         parameter.nlevel=7;
%         parameter.k_down=0.8;
%         parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;1 1 1];
%         parameter.useTop=useTop;
%         parameter.grid_space=[2 2 2;2 2 1;4 4 2;6 6 3;8 8 4; 10 10 5;12 12 6];
%         parameter.x_tot=71;
%         parameter.y_tot=71;
%         parameter.z_tot=71;
%         parameter.metric='MIND';
%         parameter.labels=[3 3 3;5 5 5;7 7 7;7 7 5;7 7 7;7 7 7;13 13 13];
%         parameter.presmooth=0;
%         parameter.spc=spc;
%         parameter.smooth_co=[1,80,10,10,1,1,0.8];
%         parameter.Tcoe_n=[20,0.001,0.001,0.001,0.05,0.05,0.005];
%         parameter.end=2;
%         parameter.top_co=[1,2.4/((12/2)^2),2.4/(((12/4)^2)*(6/2)),2.4/(((12/6)^2)*(6/3)),2.4/(((12/8)^2)*(6/4)),2.4/(((12/10)^2)*(6/5)),2.4];
%         parameter.resize=0;
%         parameter.dist_co=100;
    case 10
        
                        parameter.nlevel=6;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1;1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[2 2 2;2 2 2;4 4 4;4 4 4; 6 6 6; 8 8 8];
        parameter.x_tot=203;
        parameter.y_tot=203;
        parameter.z_tot=203;
        parameter.metric='MIND';
        parameter.labels=[3 3 3;5 5 5;5 5 5;7 7 7;7 7 7;9 9 13;];
        parameter.presmooth=0;
        parameter.spc=spc;
        if useTop==1
            parameter.smooth_co=0.05;
        else
            parameter.smooth_co=0.05;
        end
        parameter.end=6;
        parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3   ,spc(1)*spc(2)*spc(3)/(10/2)^3,  spc(1)*spc(2)*spc(3)/(10/2)^3,    spc(1)*spc(2)*spc(3)/(10/6)^3,  spc(1)*spc(2)*spc(3)/(10/8)^3,    spc(1)*spc(2)*spc(3)/1^3;];
        parameter.top_co=[1   ,0.008,  0.06,  0.125,  0.42, 1;];
        parameter.Tcoe_n=5;
        parameter.resize=1;
        parameter.dist_co=1;
        
%         parameter.nlevel=7;
%         parameter.k_down=0.8;
%         parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;1 2 2];
%         parameter.useTop=useTop;
%         parameter.grid_space=[2 2 2;2 2 1;4 4 2;6 6 3;8 8 4; 10 10 5;12 12  6];
%         parameter.x_tot=71;
%         parameter.y_tot=71;
%         parameter.z_tot=71;
%         parameter.metric='MIND';
%         parameter.labels=[3 3 3;5 5 5;7 7 7;7 7 5;7 7 7;7 7 7;13 11 11];
%         parameter.presmooth=0;
%         parameter.spc=spc;
%         parameter.smooth_co=[1,80,10,10,1,1,1];
%         parameter.Tcoe_n=[20,0.001,0.001,0.001,0.05,0.05,0.1];
%         parameter.end=2;
%         parameter.top_co=[spc(1)*spc(2)*spc(3)/1^3  ,spc(1)*spc(2)*spc(3)/1^3 ,spc(1)*spc(2)*spc(3)/(10/2)^3,   spc(1)*spc(2)*spc(3)/64,  spc(1)*spc(2)*spc(3)/(10/6)^3,   spc(1)*spc(2)*spc(3)/(10/8)^3,  spc(1)*spc(2)*spc(3)/1^3;];
%         parameter.resize=0;
%         parameter.dist_co=100;

end
end