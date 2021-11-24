clear;
load('E:\3d-feature-based\Registration-master\ExSeq\output\TPSMap_sa0916dncv_round002.mat');
r=239;c=155;s=84;
T=zeros(r,c,s,3);
for n=1:length(out1D_total)
    num=in1D_total(n);
   [x_index,y_index,z_index]=ind2sub([r,c,s],num);
   [x_trans,y_trans,z_trans]=ind2sub([r,c,s],out1D_total(n));
   T(x_index,y_index,z_index,1)=x_trans-x_index;
   T(x_index,y_index,z_index,2)=y_trans-y_index;
   T(x_index,y_index,z_index,3)=z_trans-z_index; 
end
