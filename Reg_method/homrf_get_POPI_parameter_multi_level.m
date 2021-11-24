function parameter=homrf_get_POPI_parameter_multi_level(useTop,spc)


        parameter.nlevel=5;
        parameter.k_down=0.8;
        parameter.quant=[1 1 1; 1 1 1; 1 1 1; 1 1 2; 1 1 3];
        parameter.useTop=useTop;
        parameter.grid_space=[1 1 1; 5 5 5; 4 4 4; 6 6 6; 8 8 8; ];
        parameter.x_tot=103;
        parameter.y_tot=103;
        parameter.z_tot=103;
        parameter.metric='MIND';
        parameter.labels=[5 5 5; 3 3 3; 3 3 3;5 5 5; 7 7 7;];
        parameter.presmooth=0;
        parameter.spc=spc;
        if useTop==1
            parameter.smooth_co = 0.1;
        else
            parameter.smooth_co = 0.1;
        end
        parameter.end = 3;
      
        parameter.top_co= 1;
        parameter.Tcoe_n=5;
        parameter.resize=1;
        parameter.dist_co=1;

    

end
