% This is the test of the VAMPIRE chanllenge dataset

clear;
addpath(genpath('../VAMPIRE_FullDatabase_MHA/'));

RefVI_type = 'DTPA-SPECT';
if strcmp(RefVI_type,'Galligas')
    file_path = 'H:/VAMPIRE_FullDatabase_MHA/Study01_Galligas-PET/';
    subject_num = 25;
    num_phase = 5;
else
    if strcmp(RefVI_type,'Xenon')
        file_path = 'H:/VAMPIRE_FullDatabase_MHA/Study02_Xenon-CT/';
        subject_num = 4;
        num_phase = 8;
    else
        if strcmp(RefVI_type,'DTPA-SPECT')
            file_path = 'H:/VAMPIRE_FullDatabase_MHA/Study03_DTPA-SPECT/';
            subject_num = 21;
            num_phase = 10;
        end
    end
end
for num = 17:17
    if num < 10
        CT_path = [file_path,'Subject_0',num2str(num),'/CT/'];
        Ref_path = [file_path,'Subject_0',num2str(num),'/GT/'];
    else
        CT_path = [file_path,'Subject_',num2str(num),'/CT/'];
        Ref_path = [file_path,'Subject_',num2str(num),'/GT/'];
    end
    Ref_Filename = [Ref_path,'VentImage.mha'];
    info = mha_read_header(Ref_Filename);
    vent_Ref = mha_read_volume(info);
    vol_ct = [];
    if strcmp(RefVI_type,'Xenon')
        avg_img = [];
    else
        avg_img_path = [CT_path,'AverageImage.mha'];
        info = mha_read_header(avg_img_path);
        avg_img = mha_read_volume(info);
    end
    coronal = [];
    for phase = 1:num_phase
        if phase < 10
            CT_Filename = [CT_path,'PhaseImage_0',num2str(phase),'.mha'];
        else
            CT_Filename = [CT_path,'PhaseImage_',num2str(phase),'.mha'];
        end
        info = mha_read_header(CT_Filename);
        vol_ct_cur = mha_read_volume(info);
        
        cross_section = squeeze(vol_ct_cur(:,:,20));
        coronal_cur = squeeze(vol_ct_cur(:,256,:));
        coronal_cur = imrotate(coronal_cur,90);
        coronal_cur = imresize(coronal_cur,[440,492],'bilinear');
        coronal = cat(3,coronal,coronal_cur);
        
%         subplot(1,5,phase)
% %         imshow(cross_section,[])
% %         title(['cross_',num2str(phase)]);
%         
%         imshow(coronal_cur,[])
%         title(['coronal_',num2str(phase)]);
        vol_ct = cat(4,vol_ct,vol_ct_cur);
        
    end
end
coronal_vertify = squeeze(vol_ct_cur(:,10,:));
coronal_vertify = imrotate(coronal_vertify,90);
coronal_vertify = imresize(coronal_vertify,[280,512],'bilinear');
imshowpair_xp(coronal(:,:,1),coronal(:,:,2))
% diff = mean(vol_ct,4)-avg_img;