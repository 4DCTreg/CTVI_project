function parameter=homrf_get_parameter_v1(idx,useTop,spc)
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
        parameter.nlevel=1;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1;1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[ 8 8 8];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[9 9 13;];
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

    case 3

               
        parameter.nlevel=1;
        parameter.k_down=0.8;
        parameter.quant=[ 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[ 8 8 8];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[9 9 13;];
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
    case 4

        
        parameter.nlevel=1;
        parameter.k_down=0.8;
        parameter.quant=[ 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[ 8 8 8];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[9 9 13;];
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
    case 5
        parameter.nlevel=1;
        parameter.k_down=0.8;
        parameter.quant=[ 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[ 8 8 8];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[13 13 13;];
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

    case 6

        
        parameter.nlevel=1;
        parameter.k_down=0.8;
        parameter.quant=[ 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[ 8 8 8];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[9 9 13;];
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
    case 7
        
        parameter.nlevel=1;
        parameter.k_down=0.8;
        parameter.quant=[ 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[ 6 6 6];
        parameter.x_tot=203;
        parameter.y_tot=203;
        parameter.z_tot=203;
        parameter.metric='MIND';
        parameter.labels=[9 9 13;];
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
        

    case 8

        
        
        parameter.nlevel=1;
        parameter.k_down=0.8;
        parameter.quant=[ 2 2 2];
        parameter.useTop=useTop;
        parameter.grid_space=[ 6 6 6];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[13 13 13;];
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
    case 9
        
        parameter.nlevel=1;
        parameter.k_down=0.8;
        parameter.quant=[ 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[ 6 6 6];
        parameter.x_tot=203;
        parameter.y_tot=203;
        parameter.z_tot=203;
        parameter.metric='MIND';
        parameter.labels=[9 9 13;];
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
        

    case 10
        
        parameter.nlevel=1;
        parameter.k_down=0.8;
        parameter.quant=[ 1 1 2];
        parameter.useTop=useTop;
        parameter.grid_space=[ 8 8 8];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[9 9 13;];
        parameter.presmooth=0;
        parameter.spc=spc;
        if useTop==1
            parameter.smooth_co=0.05;
        else
            parameter.smooth_co=0.05;
        end
        parameter.end=1;
      
        parameter.top_co = 1;
        parameter.Tcoe_n = 5;
        parameter.resize=1;
        parameter.dist_co=1;
        


end
end