import torch
import torch.nn.functional as F
import numpy as np
import scipy.io as scio
path = 'inter.mat'
data = scio.loadmat(path)

target_m = data['Target_r']
filter_m = data['filter']
group_m = data['group']
ngroup = int(group_m)

target_tensor = torch.from_numpy(target_m).cuda()
filter_tensor = torch.from_numpy(filter_m).cuda()
output3d = F.conv3d(target_tensor, filter_tensor, groups=ngroup)
result = output3d.cpu().numpy()
save_fname = 'results.mat'
scio.savemat(save_fname, {'array':result})




