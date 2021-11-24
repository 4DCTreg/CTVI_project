function show_supervoxel_3D(I,Label)
sup_all_sz = size(Label);
%% Show the output
Slice = round(size(Label,2)/2);
Transection = round(size(Label,3)/2);
Label_tran = squeeze(Label(:,:,Transection));
Image_tran = squeeze(I(:,:,Transection));
Image_2D = squeeze(I(:,Slice,:,1));
Image_2D = imrotate(Image_2D,90);
Label1 = squeeze(Label(:,Slice,:,1));
Label1 = imrotate(Label1,90);
k1 = unique(Label1);
k2 = unique(Label_tran);
% nor_index = [1:length(k1)];
Label_sup = zeros(size(Label1));
for nor_id = 1:length(k1)
    index_cur = find(Label1==k1(nor_id));
    Label_sup(index_cur) = nor_id;
end

Label_sup = imrotate(Label_sup,180);

Label_sup_tran = zeros(size(Label_tran));
for nor_trans = 1:length(k2)
        index_cur = find(Label_tran==k2(nor_trans));
    Label_sup_tran(index_cur) = nor_trans;
end
% Label2 = zeros(size(Image_2D));
% BW = zeros(size(Image_2D));
% BW = logical(BW);
% mycolormap = (prism(length(k1)+10));
% color_patch = zeros([size(Image_2D),3]);

% for idx = 1:numel(k1) % 1:k
%     c_k = k1(idx);
%     L = zeros(size(Image_2D));
%     L(Label1==c_k)=1;
%     BW2 = L;
%     BW_temp = edge(BW2);
%     Label2 = Label2+double(BW2)*c_k;
%     BW = BW|BW_temp;
%     temp_index = find(Label1 == c_k); 
%     [temp_x,temp_y] = ind2sub(size(Image_2D),temp_index);
%     color_index = round(rand*length(k1))+1;
%     for layer = 1:3
%         color_patch(temp_x,temp_y,1) = mycolormap(color_index,1);
%         color_patch(temp_x,temp_y,2) = mycolormap(color_index,2);
%         color_patch(temp_x,temp_y,3) = mycolormap(color_index,3);
%     end
% end

mycolormap1 = colorcube(length(k1));
mycolormap2 = colorcube(length(k2));
mycolormap = ([mycolormap1;gray(length(k1))]);
Image_tran = img_thr(Image_tran,min(Image_tran(:)),max(Image_tran(:)),1)+1;
Label_sup = img_thr(Label_sup,min(Label_sup(:)),max(Label_sup(:)),1);
% for P = 1:1
%     Image_2D = squeeze(I(:,Slice,:,P));
%     Image_2D = imrotate(Image_2D,90);
%     Image_2D_thr = img_thr(Image_2D,min(Image_2D(:)),max(Image_2D(:)),1);
%     BW_Color = repmat(Image_2D_thr,1,1,3);
%     BW_Color = uint8(BW_Color*255);
%     Image_2D_all = repmat(Image_2D_thr,1,1,3);
%     for layer = 1:2
%         tempLayer = BW_Color(:,:,layer);
%         tempLayer(BW) = 255;
%         BW_Color(:,:,layer) = tempLayer;
%     end
%     tempLayer = BW_Color(:,:,3);
%     tempLayer(BW) = 0;
%     BW_Color(:,:,3) = tempLayer;
%     figure(P);
%     subplot(1,3,1); imshow(Image_2D_all,[])
%     title(['Original CT'])
%     subplot(1,3,2); imshow(BW_Color,[])
%     title('SuperVoxel')
%     subplot(1,3,3); imagesc(Label_sup);
%     colormap(mycolormap);
%     title('Color Patch')
%     axis off
% end


k1 = Slice;
[y,z] = meshgrid(1:sup_all_sz(1),1:sup_all_sz(3));%产生网格

 x = k1*ones(size(Label_sup));
 b = surf(x,y,z,Label_sup);
%     colormap(mycolormap1);
     set(b,'linestyle','none');%隐藏网格
 hold on
 k2 = Transection;
[x,y] = meshgrid(1:sup_all_sz(2),1:sup_all_sz(1));%产生网格

 z = k2*ones(size(Label_tran,1),size(Label_tran,2));
 c = surf(x,y,z,Image_tran,'Facealpha',0.8);
 h2Data = get(c,'CData');
    colormap(mycolormap);
%     colormap(c,gray);
caxis([0,2]);
 grid off;

 set(c,'linestyle','none');%隐藏网格
%  view(-38,7);hold on;
 axis off
end