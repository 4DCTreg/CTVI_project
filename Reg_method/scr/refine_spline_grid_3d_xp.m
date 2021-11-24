function [Knots_n] = refine_spline_grid_3d_xp(Knots, new_spacing, new_vsz, O_trans_X,O_trans_Y,O_trans_Z)
ksz_new = ceil(new_vsz./new_spacing);
Kn = zeros([ksz_new, 3]);
    cur_sz = size(O_trans_X);
    for i = 1:cur_sz(1)
        for j = 1:cur_sz(2)
            for k = 1:cur_sz(3)
                i_start = O_trans_X(i,j,k,1);
                i_end = O_trans_X(i,j,k,2);
                j_start = O_trans_Y(i,j,k,1);
                j_end = O_trans_Y(i,j,k,2);
                k_start = O_trans_Z(i,j,k,1);
                k_end = O_trans_Z(i,j,k,2);
                Kn(i_start:i_end,j_start:j_end,k_start:k_end,1) = Knots(i,j,k,1);
                Kn(i_start:i_end,j_start:j_end,k_start:k_end,2) = Knots(i,j,k,2);
                Kn(i_start:i_end,j_start:j_end,k_start:k_end,3) = Knots(i,j,k,3);
                
            end
        end
    end
    
    Knots_n = Kn;

    
    

end