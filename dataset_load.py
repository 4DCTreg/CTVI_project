import os
import sys
import glob
import numpy as np

import utils

def volgen_VAMPIRE(
        vol_names,
        n_paired=6,
        batch_size=1,
        patch_size=[128, 128, 128],
        np_var='vol',
        pad_shape=None,
        resize_factor=1,
        add_feat_axis=True
):

    # convert glob path to filenames
    if isinstance(vol_names, str):
        if os.path.isdir(vol_names):
            vol_names = os.path.join(vol_names, '*')
        vol_names = glob.glob(vol_names)
    while True:
        index = np.random.randint(n_paired, size=1)
        index_ct1 = index[0] * 6
        index_ct2 = index[0] * 6 + 1
        index_ct3 = index[0] * 6 + 2
        index_ct4 = index[0] * 6 + 3
        index_ct5 = index[0] * 6 + 4
        index_gt = index[0] * 6 + 5

        load_params = dict(np_var=np_var, add_batch_axis=True, add_feat_axis=add_feat_axis, pad_shape=pad_shape,
                           resize_factor=resize_factor)
        imgs_ct1 = utils.load_volfile(vol_names[index_ct1], **load_params)
        imgs_ct2 = utils.load_volfile(vol_names[index_ct2], **load_params)
        imgs_ct3 = utils.load_volfile(vol_names[index_ct3], **load_params)
        imgs_ct4 = utils.load_volfile(vol_names[index_ct4], **load_params)
        imgs_ct5 = utils.load_volfile(vol_names[index_ct5], **load_params)
        imgs_gt = utils.load_volfile(vol_names[index_gt], **load_params)


        full_size = imgs_ct1.shape
        col_idx = np.random.randint(0, full_size[2] - patch_size[1] + 1, batch_size)
        row_idx = np.random.randint(0, full_size[1] - patch_size[0] + 1, batch_size)
        sl_idx = np.random.randint(0, full_size[3] - patch_size[2] + 1, batch_size)
        imgs_ct1_patches = [imgs_ct1[0, row_idx[i]:row_idx[i] + patch_size[0], col_idx[i]:col_idx[i] + patch_size[1], sl_idx[i]:sl_idx[i] + patch_size[2], 0].reshape(1,
                     patch_size[0], patch_size[1],  patch_size[2], 1) for i in range(batch_size)]
        imgs_ct2_patches = [imgs_ct2[0, row_idx[i]:row_idx[i] + patch_size[0], col_idx[i]:col_idx[i] + patch_size[1], sl_idx[i]:sl_idx[i] + patch_size[2], 0].reshape(1,
                     patch_size[0], patch_size[1],  patch_size[2], 1) for i in range(batch_size)]
        imgs_ct3_patches = [imgs_ct3[0, row_idx[i]:row_idx[i] + patch_size[0], col_idx[i]:col_idx[i] + patch_size[1],
                            sl_idx[i]:sl_idx[i] + patch_size[2], 0].reshape(1,
                                                                            patch_size[0], patch_size[1], patch_size[2],
                                                                            1) for i in range(batch_size)]
        imgs_ct4_patches = [imgs_ct4[0, row_idx[i]:row_idx[i] + patch_size[0], col_idx[i]:col_idx[i] + patch_size[1],
                            sl_idx[i]:sl_idx[i] + patch_size[2], 0].reshape(1,
                                                                            patch_size[0], patch_size[1], patch_size[2],
                                                                            1) for i in range(batch_size)]
        imgs_ct5_patches = [imgs_ct5[0, row_idx[i]:row_idx[i] + patch_size[0], col_idx[i]:col_idx[i] + patch_size[1],
                            sl_idx[i]:sl_idx[i] + patch_size[2], 0].reshape(1,
                                                                            patch_size[0], patch_size[1], patch_size[2],
                                                                            1) for i in range(batch_size)]
        imgs_gt_patches = [imgs_gt[0, row_idx[i]:row_idx[i] + patch_size[0], col_idx[i]:col_idx[i] + patch_size[1],
                            sl_idx[i]:sl_idx[i] + patch_size[2], 0].reshape(1,
                                                                            patch_size[0], patch_size[1], patch_size[2],
                                                                            1) for i in range(batch_size)]

        vols_input = [np.concatenate(imgs_ct1_patches, axis=0)]
        vols_input.append(np.concatenate(imgs_ct2_patches))
        vols_input.append(np.concatenate(imgs_ct3_patches))
        vols_input.append(np.concatenate(imgs_ct4_patches))
        vols_input.append(np.concatenate(imgs_ct5_patches))
        vols_gt = [np.concatenate(imgs_gt_patches, axis=0)]


        yield (vols_input, vols_gt)


def VAMPIRE_gen(vol_names, bidir=False, batch_size=1, prob_same=0, no_warp=False, **kwargs):
    """
    Generator for scan-to-scan registration.
    Parameters:
        vol_names: List of volume files to load.
        bidir: Yield input image as output for bidirectional models. Default is False.
        batch_size: Batch size. Default is 1.
        prob_same: Induced probability that source and target inputs are the same. Default is 0.
        no_warp: Excludes null warp in output list if set to True (for affine training). Default if False.
        kwargs: Forwarded to the internal volgen generator.
    """
    zeros = None
    n_paired = (len(vol_names) / 6)
    gen_paired = volgen_VAMPIRE(vol_names, int(n_paired), batch_size=batch_size, patch_size=[160, 128, 128],
                                          **kwargs)

    while True:
        scan1, scan2 = next(gen_paired)
        scan1_liver_seg = scan1[1]
        scan2_liver_seg = scan2[1]
        scan1_vessel_seg = scan1[2]
        scan2_vessel_seg = scan2[2]
        scan1_vessel_exp_seg = scan1[3]
        scan2_vessel_exp_seg = scan2[3]
        scan1_edge_seg = scan1[4]
        scan2_edge_seg = scan2[4]

        scan1 = scan1[0]
        scan2 = scan2[0]

        # some induced chance of making source and target equal
        if prob_same > 0 and np.random.rand() < prob_same:
            if np.random.rand() > 0.5:
                scan1 = scan2
            else:
                scan2 = scan1

        # cache zeros
        if not no_warp and zeros is None:
            shape = scan1.shape[1:-1]
            zeros = np.zeros((batch_size, *shape, len(shape)))

        invols = [scan1, scan2, scan1_liver_seg, scan2_liver_seg, scan1_vessel_seg, scan2_vessel_seg,
                  scan1_vessel_exp_seg, scan2_vessel_exp_seg, scan1_edge_seg, scan2_edge_seg]
        outvols = [scan2, scan1] if bidir else [scan2]
        if not no_warp:
            outvols.append(zeros)
            outvols.append(zeros)
            outvols.append(zeros)
            outvols.append(zeros)
            outvols.append(zeros)
            outvols.append(zeros)
        outvols.append(scan2_liver_seg)
        outvols.append(scan2_vessel_seg)
        outvols.append(zeros)

        yield (invols, outvols)