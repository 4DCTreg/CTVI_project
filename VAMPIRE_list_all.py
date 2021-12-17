import numpy as np
import random
random.seed(123)  # Python random module.
path_patches = 'D:/CTVI_data/Step1_pre_processing/'
train_galligas = np.random.choice(np.arange(25), size=21, replace=False)
test_galligas = np.delete(np.arange(25), train_galligas)
train_dtpa = np.random.choice(np.arange(21), size=18, replace=False)
test_dtpa = np.delete(np.arange(21), train_dtpa)
train_dtpa = train_dtpa + 25
test_dtpa = test_dtpa + 25
train_index = np.hstack((train_dtpa, train_galligas))
test_index = np.hstack((test_galligas, test_dtpa))

train_index = train_index+1
test_index = test_index+1

with open('VAMPIRE_train_all.txt', 'w') as f:
    for n_index in train_index:
        CT_1 = 'Sub_' + str(n_index) + '_CT_img_1.nii'
        CT_2 = 'Sub_' + str(n_index) + '_CT_img_2.nii'
        CT_3 = 'Sub_' + str(n_index) + '_CT_img_3.nii'
        CT_4 = 'Sub_' + str(n_index) + '_CT_img_4.nii'
        CT_5 = 'Sub_' + str(n_index) + '_CT_img_5.nii'
        GT_img = 'Sub_' + str(n_index) + '_GT_img.nii'

        f.write(path_patches + CT_1)
        f.write('\n')
        f.write(path_patches+CT_2)
        f.write('\n')
        f.write(path_patches+CT_3)
        f.write('\n')
        f.write(path_patches+CT_4)
        f.write('\n')
        f.write(path_patches+CT_5)
        f.write('\n')
        f.write(path_patches+GT_img)
        f.write('\n')

f.close()

with open('VAMPIRE_test_all.txt', 'w') as h:
    for n_index in test_index:
        CT_1 = 'Sub_' + str(n_index) + '_CT_img_1.nii'
        CT_2 = 'Sub_' + str(n_index) + '_CT_img_2.nii'
        CT_3 = 'Sub_' + str(n_index) + '_CT_img_3.nii'
        CT_4 = 'Sub_' + str(n_index) + '_CT_img_4.nii'
        CT_5 = 'Sub_' + str(n_index) + '_CT_img_5.nii'
        GT_img = 'Sub_' + str(n_index) + '_GT_img.nii'
        h.write(path_patches + CT_1)
        h.write('\n')
        h.write(path_patches + CT_2)
        h.write('\n')
        h.write(path_patches + CT_3)
        h.write('\n')
        h.write(path_patches + CT_4)
        h.write('\n')
        h.write(path_patches + CT_5)
        h.write('\n')
        h.write(path_patches + GT_img)
        h.write('\n')
h.close()