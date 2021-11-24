% XP debug
% follow the https://github.com/dgoodwin208/Registration/wiki
% First you should download the VLFeat toolbox 
clear;
addpath(genpath('E:\3d-feature-based\Registration-master\Registration\Code\'));
calculateDescriptors(1,1,9);
calculateDescriptors(2,1,9);
registerWithDescriptors_correct(2); 