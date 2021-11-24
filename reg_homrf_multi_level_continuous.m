function [Knots_n_tot,Knots_n]=reg_homrf_multi_level_continuous(parameter,...
    volmov,volfix,init_size,pts_fix_or,pts_mov_or,bszv,crop_v,Knots_n_pre)

[nlevel,k_down,x_tot, y_tot, z_tot, useTop, grid_space, metric, labels_search,...
    presmooth, quant, smooth_co_pre, endl, Tcoe_n, top_co, spc, resize, dist_co] = homrf_process_parameters_continuous_v1(parameter);

max_gs=grid_space(nlevel,:);


if presmooth == 1
    dist=Get_Smooth_tot(x_tot,y_tot,z_tot,spc);
else
    dist=[];
end

grid_sz = zeros(nlevel,3);
count_jaco = zeros(nlevel,1);
TRE_current = zeros(nlevel,1);
STD_cur = zeros(nlevel,1);
for level=nlevel:-1:endl
    
    labels=DIR_Process_labels(labels_search,level);
    
    % Current Grid Size
    cur_grid_space = grid_space(level,:);
    
    smooth_co = prod(max_gs)/prod(cur_grid_space)*smooth_co_pre;
    %     smooth_co = smooth_co_pre;
    
    % Make the registration control grid
    [O_trans_X,O_trans_Y,O_trans_Z,dx,dy,dz]=make_control_grid_3D_Odd(cur_grid_space,size(volmov));
    % Control grid dimensions (in the image)
    griddim=[size(O_trans_X,1),size(O_trans_X,2),size(O_trans_X,3)];
    grid_sz(level,:)=griddim;
    
    if level == nlevel
        Knots_n = get_Knots_pre(Knots_n_pre,O_trans_X,O_trans_Y,O_trans_Z,griddim);
        previous = get_previous_displacement(Knots_n,griddim);
    else
        Knots_n = refine_spline_grid_3D_nl(Knots_n, O_trans_X_or,...
            O_trans_Y_or,O_trans_Z_or, O_trans_X,O_trans_Y,O_trans_Z);
        previous=get_previous_displacement(Knots_n,griddim);
    end
    previous=round(previous);
    labelstot=get_labels(previous,x_tot,y_tot,z_tot);
    O_trans_X_or=O_trans_X;O_trans_Y_or=O_trans_Y;O_trans_Z_or=O_trans_Z;
    
    % Get the neighbor
    neighbor=Get_Neighbor_18(griddim(1),griddim(2),griddim(3));
    
    % Get the unary term
    tic
    unary=Get_Data_Term_3D_mask(volmov,volfix,[],[],griddim,O_trans_X,O_trans_Y,...
        O_trans_Z,labels,level,previous,quant,metric,cur_grid_space);
    toc
    
    tic
    [Ln,labelstot_Jaco] = UGM_Decode_ICM_3D(labels,dist,dist_co,neighbor,...
        griddim,unary,level,k_down,useTop,previous,x_tot,y_tot,z_tot,labelstot,...
        grid_space,spc,quant,smooth_co,Tcoe_n,top_co);
    toc
    
    Knots_min = Knots_Displacement_3D_Multi_forwards(Ln,labels.sx,labels.sy,...
        labels.sz,griddim,level,quant);
    
    Knots_n = Knots_n + Knots_min;
    
    Knots_n_tot = refine_spline_grid_3d(Knots_n,[1,1,1], size(volmov), O_trans_X, O_trans_Y, O_trans_Z);

%     if resize ==1
%         Tptv_rsz = cat(4, volresize(Knots_n_tot(:,:,:,1), bszv), ...
%             volresize(Knots_n_tot(:,:,:,2), bszv), volresize(Knots_n_tot(:,:,:,3), bszv));
%         voldef = volresize(volmov,bszv);
%     end
%     
%     [~, Knots_n_rsz] = uncrop_data(voldef, Tptv_rsz, crop_v, init_size);
    
end

end