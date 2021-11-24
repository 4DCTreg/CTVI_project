function vent_img = Ventilation_Supervoxel_v1(vol_ex,vol_in,vent_jac_interp,SLIC_Labels_3D,V_after,...
    C_p_in,labels_index,Tptv,spc)

Numlabels = size(labels_index,1);

im_sz = size(vol_ex);
ex_mask_1 = find(vol_ex>-250);
ex_mask_2 = find(vol_ex<=-1024);
ex_mask_all = [ex_mask_1;ex_mask_2];
vol_ex_sup = vol_ex;
vol_ex_sup(ex_mask_1) = -250;
vol_ex_sup(ex_mask_2) = -1024;
in_mask_1 = find(vol_in>-250);
in_mask_2 = find(vol_in<=-1024);
vol_in_sup = vol_in;
vol_in_sup(in_mask_1) = -250;
vol_in_sup(in_mask_2) = -1024;



tv_reg = 0.01;
vent_img = zeros(im_sz);
L_size = size(SLIC_Labels_3D);
for id = 1:Numlabels
    vol_ex_index = find(SLIC_Labels_3D==labels_index(id));
    
    sup_area_mask = zeros(L_size(1),L_size(2),L_size(3));
    sup_area_mask(vol_ex_index) = 1;
    
    [c_x,c_y,c_z] = ind2sub(L_size,vol_ex_index);
    
    
    
    new_x = c_x ;
    new_y = c_y ;
    new_z = c_z ;
    
    
    %     figure;
    %     plot(shp);
    min_x = floor(min(new_x));
    if min_x <= 0
        min_x = 1;
    end
    max_x = ceil(max(new_x));
    if max_x > L_size(1)
        max_x = L_size(1);
    end
    
    min_y = floor(min(new_y));
    if min_y <= 0
        min_y = 1;
    end
    
    max_y = ceil(max(new_y));
    if max_y > L_size(2)
        max_y = L_size(2);
    end
    min_z = floor(min(new_z));
    if min_z<=0
        min_z = 1;
    end
    max_z = ceil(max(new_z));
    if max_z > L_size(3)
        max_z = L_size(3);
    end
    [Y_grid,X_grid,Z_grid] = meshgrid(min_y:max_y,min_x:max_x,min_z:max_z);
    
    X_index_cur = reshape(X_grid,numel(X_grid),1);
    Y_index_cur = reshape(Y_grid,numel(Y_grid),1);
    Z_index_cur = reshape(Z_grid,numel(Z_grid),1);
    area_index = sub2ind(L_size,X_index_cur,Y_index_cur,Z_index_cur);
    if size(unique(c_x),1)==1||size(unique(c_y),1)==1||size(unique(c_z),1)==1
        C_p_ex_cur = [];
        supvoxel_sz_ex = zeros(1,6);
    else
        C_p_ex_cur = sup_area_mask(min_x:max_x,min_y:max_y,min_z:max_z);
        C_p_ex_cur = reshape(C_p_ex_cur,numel(X_grid),1);
        C_p_ex_tv = cat(2,area_index,C_p_ex_cur);
        supvoxel_sz_ex =[min_x,max_x,min_y,max_y,min_z,max_z];
    end
    
    volin_sup = C_p_in{id};
    if isempty(volin_sup)||isempty(C_p_ex_cur)
        continue
    else
        volin_index = volin_sup(:,1);
        volin_index_mask = volin_sup(:,2);
        volex_air = (-vol_ex_sup(vol_ex_index));
        volex_air_sum = sum(volex_air(:));
        
        volin_air = ((-vol_in_sup(volin_index))).*volin_index_mask;
        volin_air_sum = sum(volin_air(:));
%         sup_vent = (volin_air_sum - volex_air_sum)/sum(volex_air_sum)*size(vol_ex_index,1)^2;
        sup_vent = (volin_air_sum - volex_air_sum);
        

        
        n_voxel = size(find(C_p_ex_cur==1),1);
        vent_voxel = zeros(n_voxel,1);
        sup_jac_coef = vent_jac_interp(area_index);
        sup_no_mask = find(C_p_ex_cur==0);
        sup_jac_coef(sup_no_mask)=[];
        sup_coef_min = min(sup_jac_coef);
        sup_coef_max = max(sup_jac_coef);
        if ceil(sup_coef_min)==ceil(sup_coef_max)
            sup_jac_coef = ones(size(sup_jac_coef,1),1);
        else
            
            sup_jac_coef = 1./sup_jac_coef;
            sup_jac_coef = img_thr_v1(sup_jac_coef, sup_coef_min, sup_coef_max, 5);
            sup_vent_index = sup_jac_coef;
            sup_vent_index_sort = sort(sup_vent_index);
            min_v = sup_vent_index_sort(ceil(0.02*numel(sup_jac_coef)));
            max_v = sup_vent_index_sort(floor(0.98*numel(sup_jac_coef)));
            
            %             max_v = max(sup_jac_coef);
            %             min_v = min(sup_jac_coef);
            if ceil(max_v) == ceil(min_v)
                sup_jac_coef = ones(size(sup_jac_coef,1),1);
            else
                                sup_jac_coef = img_thr_v1(sup_jac_coef, min(sup_jac_coef), max(sup_jac_coef), 5);
%                 sup_jac_coef = ones(size(sup_jac_coef,1),1);
            end
        end
        
        % sup_jac_coef = sup_jac_coef/sum(sup_jac_coef);
        vol_sup_ex_index = area_index;
        vol_sup_ex_index(sup_no_mask) = [];
        %     sup_coef_min = min(sup_jac_coef);
        %     sup_coef_max = max(sup_jac_coef);
        %     sup_jac_coef = img_thr(sup_jac_coef, sup_coef_min, sup_coef_max, 5);
        sp_sz = supvoxel_sz_ex;
        objf = @(X)tv_grad(X, sup_vent, sup_jac_coef,C_p_ex_tv,sp_sz,tv_reg);
        options = [];
        options.display = 'off';
        % options.LS_init = 4;
        options.MaxIter = 80;
        options.MaxFunEvals = round(80 * 1.3);
        options.Method = 'lbfgs';
        
        options.DerivativeCheck = 'off';
        options.Damped = 0;
        options.LS_type = 0; options.LS_init = 8; %
        %     options.LS_type = 0; options.LS_init = 2; %
        %     options.LS_type = 1; options.LS_init = 2; %
        options.optTol = 1e-8;
        options.progTol = 1e-4;
        
        options.Corr = 85;
        if numel(vent_voxel) > (100^3)
            options.Corr = 30;
            %         options.Corr
        end
        %     options.progTol = 1e-30;
        %     options.optTol = 1e-30;
        
        options.A_step_length = 1/10;
        [xmin, fmin, ~, outp] = minFuncJoker(objf, vent_voxel(:), options);
        if sum(xmin) == 0
            ncoef = sup_vent;
        else
            ncoef = sup_vent/sum(xmin);
        end
        xmin = xmin*ncoef;
        vent_img(vol_sup_ex_index) = sup_vent;
    end
end

vent_img_mask = vent_img;
% vent_index = reshape(vent_img_mask,[numel(vol_ex),1]);
% vent_index_sort = sort(vent_index);
% min_p_value = vent_index_sort(round(0.05*numel(vol_ex)));
% max_p_value = vent_index_sort(round(0.95*numel(vol_ex)));
% vent_img_mask(find(vent_img_mask < min_p_value)) = min_p_value;
% vent_img_mask(find(vent_img_mask > max_p_value)) = max_p_value;

min_p_value = min(vent_img_mask(:));
max_p_value = max(vent_img_mask(:));
vent_img = img_thr(vent_img_mask, min_p_value, max_p_value, 1500);
vent_img(ex_mask_all)=0;




    function [xmin, fmin, v3, outp] = minFuncJoker(objf, x0, options)
        Method = getoptions(options, 'Method');
        if strcmp(Method, 'lbfgs+cg')
            options.Method = 'lbfgs';
            [xmin, fmin, v3, outp] = minFunc(objf, x0(:), options);
            options.Method = 'cg';
            [xmin, fmin, v3, outp] = minFunc(objf, xmin(:), options);
        elseif strcmp(Method, 'lbfgs+csd')
            options.Method = 'lbfgs';
            [xmin, fmin, v3, outp] = minFunc(objf, x0(:), options);
            options.Method = 'csd';
            [xmin, fmin, v3, outp] = minFunc(objf, xmin(:), options);
        elseif strcmp(Method, 'adam')
            v3=0;
            outp = 0;
            sOpt = optimset('fmin_adam');
            sOpt.GradObj = 'on';
            sOpt.MaxFunEvals = options.MaxFunEvals;
            sOpt.MaxIter = options.MaxIter;
            sOpt.Display = 'off';
            step_length = options.A_step_length;
            [xmin, fmin] = fmin_adam(objf, x0(:), step_length, [], [], [], 1, sOpt);
        else
            [xmin, fmin, v3, outp] = minFunc(objf, x0(:), options);
        end
    end
end