
clear;
addpath(genpath('isoPTVpackage/'));

basepth = 'E:/data_prj/dir_dataset/DIR_files/';

for idx = 4:4
    
    pts_struct = DIR_get_all_points_for_the_case(idx, basepth);
    [volmov, spc] = read_DIR_volume_4dCT(idx, 0, basepth);
    volmov = double(volmov);
    T00_coronal = squeeze(volmov(:,128,:));
    T00_coronal = imrotate(T00_coronal,270);
    T00_coronal = imresize(T00_coronal,[248,289],'bilinear');
    subplot(1,3,1)
    imshow(T00_coronal,[]);
    title(['T00']);
    pts_mov = pts_struct.extreme.b;
    [volfix, spc] = read_DIR_volume_4dCT(idx, 5, basepth);
    volfix = double(volfix);
    pts_fix = pts_struct.extreme.e;
    T50_coronal = squeeze(volfix(:,128,:));
    T50_coronal = imrotate(T50_coronal,270);
    T50_coronal = imresize(T50_coronal,[248,289],'bilinear');
    subplot(1,3,2)
    imshow(T50_coronal,[]);
    title(['T50'])
    subplot(1,3,3)
    imshowpair_xp(T00_coronal,T50_coronal);

end