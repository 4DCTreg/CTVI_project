load('homrf_nan.mat');
Numlabels = size(labels_index,1);
im_sz = size(vol_ex);
ex_mask_1 = find(vol_ex>0);
ex_mask_2 = find(vol_ex<=-1000);
ex_mask_all = [ex_mask_1;ex_mask_2];
vol_ex_sup = vol_ex;
vol_ex_sup(ex_mask_1) = 0;
vol_ex_sup(ex_mask_2) = -1000;
in_mask_1 = find(vol_in>0);
in_mask_2 = find(vol_in<=-1000);
vol_in_sup = vol_in;
vol_in_sup(in_mask_1) = 0;
vol_in_sup(in_mask_2) = -1000;

tv_reg = 0.11;
vent_img = zeros(im_sz);

for id = 1:Numlabels
    vol_ex_index = find(SLIC_Labels_3D==labels_index(id));
    volin_sup = C_p_in{id};
    C_p_ex_cur = C_p_ex{id};
    if isempty(volin_sup)||isempty(C_p_ex_cur) 
        continue
    else
        volin_index = volin_sup(:,1);
        volin_index_mask = volin_sup(:,2);
        volex_air = (1000-vol_ex_sup(vol_ex_index));
        volex_air_sum = sum(volex_air(:));
        volin_air = ((1000-vol_in_sup(volin_index))).*volin_index_mask;
        volin_air_sum = sum(volin_air(:));
        sup_vent = volin_air_sum - volex_air_sum;
        
        n_voxel = size(find(C_p_ex_cur(:,2)==1),1);
        vent_voxel = zeros(n_voxel,1);
        sup_jac_coef = vent_jac_interp(C_p_ex_cur(:,1));
        sup_no_mask = find(C_p_ex_cur(:,2)==0);
        sup_jac_coef(sup_no_mask)=[];
        sup_coef_min = min(sup_jac_coef);
        sup_coef_max = max(sup_jac_coef);
        if sup_coef_min==sup_coef_max
            sup_jac_coef = 2.5*ones(size(sup_jac_coef,1),1);
        else
        sup_jac_coef = img_thr(sup_jac_coef, sup_coef_min, sup_coef_max, 5);
        end
        sup_jac_coef = sup_jac_coef+1;
        vol_sup_ex_index = C_p_ex_cur(:,1);
        vol_sup_ex_index(sup_no_mask) = [];
        %     sup_coef_min = min(sup_jac_coef);
        %     sup_coef_max = max(sup_jac_coef);
        %     sup_jac_coef = img_thr(sup_jac_coef, sup_coef_min, sup_coef_max, 5);
        sp_sz = supvoxel_sz_ex(id,:);
        objf = @(X)tv_grad(X, sup_vent, sup_jac_coef,C_p_ex_cur,sp_sz,tv_reg);
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
        
        options.A_step_length = 1/2.4;
        [xmin, fmin, ~, outp] = minFuncJoker(objf, vent_voxel(:), options);
        if sum(xmin)==0
            ncoef = 1;
        else
            ncoef = sup_vent/sum(xmin);
        end
        xmin = xmin*ncoef;
        
        vent_img(vol_sup_ex_index) = xmin;
    end
end
vent_img_mask = vent_img;
% vent_index = reshape(vent_img_mask,[numel(vol_ex),1]);
% vent_index_sort = sort(vent_index);
% min_p_value = vent_index_sort(round(0.01*numel(vol_ex)));
% max_p_value = vent_index_sort(round(0.99*numel(vol_ex)));
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
