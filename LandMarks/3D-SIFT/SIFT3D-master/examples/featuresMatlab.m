% featuresMatlab.m
%
% This script shows how to extract SIFT3D features from a volumetric image.
%
% Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
clear;
addpath(genpath('G:\SIFT3D\SIFT3D 1.4.5\'));
% Load the image
[im, units] = imRead3D('data/1.nii.gz');
units = [0.97,0.97,2];
% Detect keypoints
keys = detectSift3D(im, 'units', units);

% Extract descriptors
[desc, coords] = extractSift3D(keys);

% Clear MEX memory
clear mex