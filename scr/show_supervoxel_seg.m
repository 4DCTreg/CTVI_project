function show_supervoxel_seg(I,Label)
%% Show the output
Slice = round(size(Label,2)/2);
Image_2D = squeeze(I(:,Slice,:,1));
Image_2D = imrotate(Image_2D,90);
Label1 = squeeze(Label(:,Slice,:,1));
Label1 = imrotate(Label1,90);
k1 = unique(Label1);
% nor_index = [1:length(k1)];
Label_sup = zeros(size(Label1));
for nor_id = 1:length(k1)
    index_cur = find(Label1==k1(nor_id));
    Label_sup(index_cur) = nor_id;
end
Label2 = zeros(size(Image_2D));
BW = zeros(size(Image_2D));
BW = logical(BW);
mycolormap = (prism(length(k1)+10));
color_patch = zeros([size(Image_2D),3]);

for idx = 1:numel(k1) % 1:k
    c_k = k1(idx);
    L = zeros(size(Image_2D));
    L(Label1==c_k)=1;
    BW2 = L;
    BW_temp = edge(BW2);
    Label2 = Label2+double(BW2)*c_k;
    BW = BW|BW_temp;
    temp_index = find(Label1 == c_k);
    [temp_x,temp_y] = ind2sub(size(Image_2D),temp_index);
    color_index = round(rand*length(k1))+1;
    for layer = 1:3
        color_patch(temp_x,temp_y,1) = mycolormap(color_index,1);
        color_patch(temp_x,temp_y,2) = mycolormap(color_index,2);
        color_patch(temp_x,temp_y,3) = mycolormap(color_index,3);
    end
end

mycolormap = colorcube(length(k1)+20);
for P = 1:1
    Image_2D = squeeze(I(:,Slice,:,P));
    Image_2D = imrotate(Image_2D,90);
    Image_2D_thr = img_thr(Image_2D,min(Image_2D(:)),max(Image_2D(:)),1);
    BW_Color = repmat(Image_2D_thr,1,1,3);
    BW_Color = uint8(BW_Color*255);
    Image_2D_all = repmat(Image_2D_thr,1,1,3);
    
    tempLayer = BW_Color(:,:,1);
    tempLayer(BW) = 255;
    BW_Color(:,:,1) = tempLayer;
    
    tempLayer = BW_Color(:,:,2);
    tempLayer(BW) = 46;
    BW_Color(:,:,2) = tempLayer;
    
    tempLayer = BW_Color(:,:,3);
    tempLayer(BW) = 46;
    BW_Color(:,:,3) = tempLayer;
    
    figure(P);
    subplot(1,3,1); imshow(Image_2D_all,[])
    title(['Original CT'])
    subplot(1,3,2); imshow(BW_Color,[])
    title('SuperVoxel')
    subplot(1,3,3); imagesc(Label_sup);
    colormap(mycolormap);
    title('Color Patch')
    axis off
end
end