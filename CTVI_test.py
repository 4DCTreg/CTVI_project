# -*- coding:utf-8 -*-
# Created Time: 2018/05/11 11:50:23
# Author: Taihong Xiao <xiaotaihong@126.com>

from dataset_load import volgen_VAMPIRE, volgen_VAMPIRE_val
from nets import Generator, Discriminator, Generator_Unet
from dataset import config
import os, argparse
import torch
import numpy as np
import scipy.io as sio
from tensorboardX import SummaryWriter
from itertools import chain
import utils
import re


class _3DGAN(object):
    def __init__(self, args, test_files, config=config):
        self.args = args
        self.attribute = args.attribute
        self.gpu = args.gpu
        self.mode = args.mode
        self.restore = args.restore

        # init dataset and networks
        self.config = config
        #self.dataset = ShapeNet(self.attribute)
        self.dataset = volgen_VAMPIRE(test_files)
        self.val_dataset = volgen_VAMPIRE_val(test_files)
        self.test_files = test_files
        self.G = Generator_Unet(5, 1)
        #self.G = Generator()
        self.D = Discriminator()

        self.adv_criterion = torch.nn.BCELoss()

        self.set_mode_and_gpu()
        self.restore_from_file()

    def set_mode_and_gpu(self):
        if self.mode == 'train':
            self.G.train()
            self.D.train()
            if self.gpu:
                with torch.cuda.device(self.gpu[0]):
                    self.G.cuda()
                    self.D.cuda()
                    self.adv_criterion.cuda()

            if len(self.gpu) > 1:
                self.G = torch.nn.DataParallel(self.G, device_ids=self.gpu)
                self.D = torch.nn.DataParallel(self.D, device_ids=self.gpu)

        elif self.mode == 'test':
            self.G.eval()
            self.D.eval()
            if self.gpu:
                with torch.cuda.device(self.gpu[0]):
                    self.G.cuda()
                    self.D.cuda()

            if len(self.gpu) > 1:
                self.G = torch.nn.DataParallel(self.G, device_ids=self.gpu)
                self.D = torch.nn.DataParallel(self.D, device_ids=self.gpu)


    def restore_from_file(self):
        if self.restore is not None:
            ckpt_file_G = os.path.join(self.config.model_dir, 'G_iter_{:06d}.pth'.format(self.restore))
            assert os.path.exists(ckpt_file_G)
            self.G.load_state_dict(torch.load(ckpt_file_G))

            if self.mode == 'train':
                ckpt_file_D = os.path.join(self.config.model_dir, 'D_iter_{:06d}.pth'.format(self.restore))
                assert os.path.exists(ckpt_file_D)
                self.D.load_state_dict(torch.load(ckpt_file_D))

            self.start_step = self.restore + 1
        else:
            self.start_step = 1

    def save_log(self):
        scalar_info = {
            'loss_D': self.loss_D,
            'loss_G': self.loss_G,
            'G_lr'  : self.G_lr_scheduler.get_last_lr()[0],
            'D_lr'  : self.D_lr_scheduler.get_last_lr()[0],
        }
        for key, value in self.G_loss.items():
            scalar_info['G_loss/' + key] = value

        for key, value in self.D_loss.items():
            scalar_info['D_loss/' + key] = value

        for tag, value in scalar_info.items():
            self.writer.add_scalar(tag, value, self.step)

    def save_img(self, save_num=5):
        for i in range(save_num):
            mdict = {
                'instance': self.fake_X[i,0].data.cpu().numpy()
            }
            sio.savemat(os.path.join(self.config.img_dir, '{:06d}_{:02d}.mat'.format(self.step, i)), mdict)

    def save_model(self):
        torch.save({key: val.cpu() for key, val in self.G.state_dict().items()}, os.path.join(self.config.model_dir, 'G_iter_{:06d}.pth'.format(self.step)))
        torch.save({key: val.cpu() for key, val in self.D.state_dict().items()}, os.path.join(self.config.model_dir, 'D_iter_{:06d}.pth'.format(self.step)))

    def train(self):
        self.writer = SummaryWriter(self.config.log_dir)
        self.opt_G = torch.optim.Adam(self.G.parameters(), lr=self.config.G_lr, betas=(0.5, 0.999))
        self.opt_D = torch.optim.Adam(self.D.parameters(), lr=self.config.D_lr, betas=(0.5, 0.999))
        self.G_lr_scheduler = torch.optim.lr_scheduler.StepLR(self.opt_G, step_size=self.config.step_size, gamma=self.config.gamma)
        self.D_lr_scheduler = torch.optim.lr_scheduler.StepLR(self.opt_D, step_size=self.config.step_size, gamma=self.config.gamma)

        initial_epoch = 0
        epochs = 200
        # start training
        for epoch in range(initial_epoch, epochs):

            epoch_total_loss_D = []
            epoch_total_loss_G = []


            for step in range(self.start_step, 1 + self.config.max_iter):
                self.step = step
                self.G_lr_scheduler.step()
                self.D_lr_scheduler.step()

                #self.real_X = next(self.dataset.gen(True))
                inputs, real_X = next(self.dataset)
                inputs = [torch.from_numpy(d).to('cuda').float().permute(0, 4, 1, 2, 3) for d in inputs]
                self.real_X = [torch.from_numpy(d).to('cuda').float().permute(0, 4, 1, 2, 3) for d in real_X]

                # self.noise = torch.randn(self.config.nchw[0], 200)
                # if len(self.gpu):
                #     with torch.cuda.device(self.gpu[0]):
                #         self.real_X = self.real_X.cuda()
                #         self.noise  = self.noise.cuda()

                self.fake_X = self.G(*inputs)

                # update D
                self.D_real = self.D(self.real_X[0])
                self.D_fake = self.D(self.fake_X.detach())
                self.D_loss = {
                    'adv_real': self.adv_criterion(self.D_real, torch.ones_like(self.D_real)),
                    'adv_fake': self.adv_criterion(self.D_fake, torch.zeros_like(self.D_fake)),
                }
                self.loss_D = sum(self.D_loss.values())

                self.opt_D.zero_grad()
                self.loss_D.backward()
                self.opt_D.step()

                # update G
                self.D_fake = self.D(self.fake_X)
                self.G_loss = {
                    'adv_fake': self.adv_criterion(self.D_fake, torch.ones_like(self.D_fake))
                }
                self.loss_G = sum(self.G_loss.values())
                self.opt_G.zero_grad()
                self.loss_G.backward()
                self.opt_G.step()

                epoch_total_loss_D.append(self.loss_D.data.cpu().numpy().item())
                epoch_total_loss_G.append(self.loss_G.data.cpu().numpy().item())

            print('epoch: {:06d}, loss_D: {:.6f}, loss_G: {:.6f}'.format(epoch, np.mean(epoch_total_loss_G), np.mean(epoch_total_loss_G)))

            if self.step % 100 == 0:
                self.save_log()

            if self.step % 1000 == 0:
                #self.save_img()
                self.save_model()

        print('Finished training!')
        self.writer.close()


    def test(self):

        self.opt_G = torch.optim.Adam(self.G.parameters(), lr=self.config.G_lr, betas=(0.5, 0.999))
        self.opt_D = torch.optim.Adam(self.D.parameters(), lr=self.config.D_lr, betas=(0.5, 0.999))
        self.G_lr_scheduler = torch.optim.lr_scheduler.StepLR(self.opt_G, step_size=self.config.step_size, gamma=self.config.gamma)
        self.D_lr_scheduler = torch.optim.lr_scheduler.StepLR(self.opt_D, step_size=self.config.step_size, gamma=self.config.gamma)
        save_path = 'D:/CTVI_data/Gan_results/'
        # fixed, fixed_affine = utils.load_volfile('D:/CTVI_data/Step1_pre_processing/Sub_1_GT_img.nii', add_batch_axis=True, ret_affine=True)
        for n_index in range(7):
            cur_name = self.test_files[n_index*6+5]
            rs = re.sub('^D:/CTVI_data/Step1_pre_processing/', '', cur_name)
            cur_save_name = save_path + rs

            self.G_lr_scheduler.step()
            self.D_lr_scheduler.step()

            #self.real_X = next(self.dataset.gen(True))
            inputs, real_X = next(self.val_dataset)
            inputs = [torch.from_numpy(d).to('cuda').float().permute(0, 4, 1, 2, 3) for d in inputs]
            self.real_X = [torch.from_numpy(d).to('cuda').float().permute(0, 4, 1, 2, 3) for d in real_X]

            # self.noise = torch.randn(self.config.nchw[0], 200)
            # if len(self.gpu):
            #     with torch.cuda.device(self.gpu[0]):
            #         self.real_X = self.real_X.cuda()
            #         self.noise  = self.noise.cuda()
            with torch.no_grad():
                self.fake_X = self.G(*inputs)

            img_g = self.fake_X
            img_g = img_g.detach().cpu().numpy().squeeze()

            utils.save_volfile_VAMPIRE(img_g, cur_save_name)

            # # update D
            # self.D_real = self.D(self.real_X[0])
            # self.D_fake = self.D(self.fake_X.detach())
            # self.D_loss = {
            #     'adv_real': self.adv_criterion(self.D_real, torch.ones_like(self.D_real)),
            #     'adv_fake': self.adv_criterion(self.D_fake, torch.zeros_like(self.D_fake)),
            # }
            # self.loss_D = sum(self.D_loss.values())
            #
            # self.opt_D.zero_grad()
            # self.loss_D.backward()
            # self.opt_D.step()
            #
            # # update G
            # self.D_fake = self.D(self.fake_X)
            # self.G_loss = {
            #     'adv_fake': self.adv_criterion(self.D_fake, torch.ones_like(self.D_fake))
            # }
            # self.loss_G = sum(self.G_loss.values())
            # self.opt_G.zero_grad()
            # self.loss_G.backward()
            # self.opt_G.step()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--attribute', type=str, help='Specify category for training.')
    parser.add_argument('-g', '--gpu', default=[0], nargs='+', type=int, help='Specify GPU ids.')
    parser.add_argument('-r', '--restore', default='000090', action='store', type=int, help='Specify checkpoint id to restore.')
    parser.add_argument('-m', '--mode', default='test', type=str, choices=['train', 'test'])
    parser.add_argument('--img-list', default='VAMPIRE_test_all.txt', help='line-seperated list of training files')
    args = parser.parse_args()
    print(args)
    test_files = utils.read_file_list(args.img_list)

    model = _3DGAN(args, test_files)
    if args.mode == 'test':
        model.test()

