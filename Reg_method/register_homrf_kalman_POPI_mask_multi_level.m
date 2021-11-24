function [Knots_n_tot,Knots_n]=register_homrf_kalman_POPI_mask_multi_level(parameter,...
    volmov,volfix,mov_mask,fix_mask,init_size,pts_fix_or,pts_mov_or,Knots_n_pre,bszv,crop_v)


nlevel=parameter.nlevel;
k_down=parameter.k_down;
x_tot=parameter.x_tot;
y_tot=parameter.y_tot;
z_tot=parameter.z_tot;
useTop=parameter.useTop;
grid_space=parameter.grid_space;
metric=parameter.metric;
labels_search=parameter.labels;
presmooth=parameter.presmooth;
quant=parameter.quant;
smooth_co_pre=parameter.smooth_co;
endl=parameter.end;
Tcoe_n=parameter.Tcoe_n;
top_co=parameter.top_co;
spc_or = parameter.spc;
% spc_orig = spc;
spc = [1,1,1];
resize=parameter.resize;
dist_co=parameter.dist_co;
max_gs=grid_space(nlevel,:);
% volsz = size(volmov);
ir = init_size(1);
ic = init_size(2);
is = init_size(3);
[Xgrid,Ygrid,Zgrid] = ndgrid(1:ir,1:ic,1:is);
Xgrid = Xgrid.*spc_or(1);
Ygrid = Ygrid.*spc_or(2);
Zgrid = Zgrid.*spc_or(3);

if presmooth==1
    dist=Get_Smooth_tot(x_tot,y_tot,z_tot,spc);     
else
    dist=[];
end

    grid_sz=zeros(nlevel,3);
    count_jaco=zeros(nlevel,1);
    TRE_current=zeros(nlevel,1);
    STD_cur=zeros(nlevel,1);
for level=nlevel:-1:endl
    
    labels=DIR_Process_labels(labels_search,level);
        
    % Current Grid Size
    cur_grid_space = grid_space(level,:);
    
    smooth_co = prod(max_gs)/prod(cur_grid_space)*smooth_co_pre;
%     smooth_co = smooth_co_pre;
    

    % Make the registration control grid 
%     [O_trans_X,O_trans_Y,O_trans_Z,dx,dy,dz]=make_control_grid_3D_homrf(cur_grid_space,size(volmov));
    [O_trans_X,O_trans_Y,O_trans_Z,dx,dy,dz]=make_control_grid_3D_Odd(cur_grid_space,size(volmov));
    % Control grid dimensions (in the image)
    griddim=[size(O_trans_X,1),size(O_trans_X,2),size(O_trans_X,3)];
    grid_sz(level,:)=griddim;

    if level==nlevel
%         Knots_n=Knots_n_pre;
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
%      neighbor=Get_Neighbor_6(griddim(1),griddim(2),griddim(3));
    % Get the unary term
    tic
    unary=Get_Data_Term_3D_mask(volmov,volfix,mov_mask,fix_mask,griddim,O_trans_X,O_trans_Y,...
        O_trans_Z,labels,level,previous,quant,metric,cur_grid_space);
    toc
    
    tic
    [Ln,labelstot_Jaco] = UGM_Decode_ICM_3D(labels,dist,dist_co,neighbor,...
        griddim,unary,level,k_down,useTop,previous,x_tot,y_tot,z_tot,labelstot,...
        grid_space,spc,quant,smooth_co,Tcoe_n,top_co);
    toc

    Knots_min=Knots_Displacement_3D_Multi_forwards(Ln,labels.sx,labels.sy,...
        labels.sz,griddim,level,quant);        
    Knots_n = Knots_n+Knots_min;

    Knots_n_tot = refine_spline_grid_3d(Knots_n,[1,1,1], size(volmov), O_trans_X,O_trans_Y,O_trans_Z);
%     Knots_n_tot = refine_spline_grid_3d_xp(Knots_n,[1,1,1], size(volmov), O_trans_X,O_trans_Y,O_trans_Z);
    if resize ==1
            Tptv_rsz = cat(4, volresize(Knots_n_tot(:,:,:,1), bszv), ...
       volresize(Knots_n_tot(:,:,:,2), bszv), volresize(Knots_n_tot(:,:,:,3), bszv));  
            voldef = volresize(volmov,bszv);
    end
        [~, Knots_n_rsz] = uncrop_data(voldef, Tptv_rsz, crop_v, init_size);

        HOMRFgriddedInterpolantX = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_rsz(:,:,:,1));
        HOMRFgriddedInterpolantY = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_rsz(:,:,:,2));
        HOMRFgriddedInterpolantZ = griddedInterpolant(Xgrid,Ygrid,Zgrid,Knots_n_rsz(:,:,:,3));
        displhomrf(:,1) = HOMRFgriddedInterpolantX(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displhomrf(:,2) = HOMRFgriddedInterpolantY(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        displhomrf(:,3) = HOMRFgriddedInterpolantZ(pts_fix_or(:,1),pts_fix_or(:,2),pts_fix_or(:,3));
        warpedLMsHOMRF = displhomrf+pts_fix_or;

        tre = mean(sqrt(sum(((warpedLMsHOMRF-pts_mov_or)).^2,2)));
        TRE_current(level) = tre;



end

end