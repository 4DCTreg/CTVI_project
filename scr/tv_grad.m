function [err, grad, fData, fReg] = tv_grad(X,sup_vent,sup_jac_coef,C_p_ex, sp_sz,tv_reg)
err = 0;
fData = 0;
fReg = 0;
sz = [sp_sz(2)-sp_sz(1)+1,sp_sz(4)-sp_sz(3)+1,sp_sz(6)-sp_sz(5)+1];
sp_cube = zeros(sz);
mask = C_p_ex(:,2);
index_mask = find(mask==1);
n_X = size(index_mask,1);
sp_cube(index_mask) = X;
diff = sum(X.*sup_jac_coef) - sup_vent;
err_Data = abs(diff);
fData = fData+err_Data;

err = fData;
g_data = sign(diff);
df = abs(sup_jac_coef)*g_data;
[f_iso, g_iso] = isotropic_tv_reg(sp_cube, [1,1,1], [], 1, 5e-3);

fReg = f_iso;
gG = g_iso(index_mask);
err = err + tv_reg*fReg;
grad = df + tv_reg*gG;