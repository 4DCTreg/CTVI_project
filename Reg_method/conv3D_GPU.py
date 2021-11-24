import torch
import torch.nn.functional as F
import numpy as np


def conv3D_GPU(img, filter, groups,channel,Nc,Nr,Ns,Nfc,Nfr,Nfs):
    img_1 = np.array(img)
    filter_1 = np.array(filter)
    groups_1 = np.array(groups)
    return filter_1





