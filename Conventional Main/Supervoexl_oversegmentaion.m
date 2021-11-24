% SuperVoxel multi-modal segmentation

% Introduction
% The algorithm used in this code is the implementation of our paper [1]
% which was fulfilled in the Lab of Vision Engineering at the University of
% Lincoln (UK). It is a modification of the method Simple Linear Iterative
% Clustering (SLIC) which was proposed by Achanta et al. [2].
% Our method is optimized for medical images such as MRI, CT, etc. The
% contributions of our codes compared to conventional 2D and 3D superpixel
% are as follows:
% ?Multi-modal input (works for single-modal, as well)
% ?Taking the spatial resolution of the medical images into account, i.e.
% the voxel resolution in X and Y directions and the slice thickness.

%% Initialization
function SLIC_Labels_3D = Supervoexl_oversegmentaion(vol_ex,spc,n_seeds)
im_sz = size(vol_ex);
im_sz_phy = im_sz.*spc;
path_phy = nthroot(((im_sz_phy(1) * im_sz_phy(2) * im_sz_phy(3))/n_seeds),3);
voxel_X = floor(path_phy);
voxel_Y = floor(path_phy);
voxel_Z = floor(path_phy);
% Load the parametrs ofrm th file "Opts.m"
Iterations = 10;
Compactness = 0.2;  % Compactness of the superpixel is equivalent to m
% voxel_X = 10;
% voxel_Y = 10;
% voxel_Z = 4;

VoxelSize = spc;
term_x = 1/VoxelSize(1);
term_y = 1/VoxelSize(2);
term_z = 1/VoxelSize(3);
Integer = floor(Compactness);
Fraction = Compactness-Integer;
Cmpt = [num2str(Integer),'_',num2str(Fraction)];

% Data and results addresses
% Case = '01';
% ProtocolList = {'FLAIR'}; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InputImage_Path = [cd,'\Data\SG_Normalized'];
% Output_Path = [cd,'\Results\SLIC'];
% mkdir(Output_Path)
% 
% for P = 1:numel(ProtocolList)
%     
%     % Selecc the protocol
%     Protocol = ProtocolList{P};
%     
%     % Displaying the patient case and protocols
%     disp(['Case:  ',Case,'  Protocol: ',Protocol])
%     
%     % Load the data
%     ImageFileName = ['SGUL_Normalized_Case_',Case,'_',Protocol,'.mat'];
%     load(fullfile(InputImage_Path,ImageFileName));
%     Image3D = Image3D_Normalized;
%     Image3D_temp = Image3D;
%     Min_I = min(Image3D_temp(:));
%     Max_I = max(Image3D_temp(:));
%     Image3D_temp = (Image3D_temp-Min_I)/(Max_I-Min_I);
%     I(:,:,:,P) = Image3D_temp;
% end
I = vol_ex;
% Initial arrays
I_temp = I(:,:,:,1);
S = nthroot((voxel_X * voxel_Y * voxel_Z),3);
X_window = floor(size(I_temp,1)/voxel_X);
Y_window = floor(size(I_temp,2)/voxel_Y);
Z_window = floor(size(I_temp,3)/voxel_Z);
r1 = ceil(voxel_X/2);
r2 = ceil(voxel_Y/2);
r3 = ceil(voxel_Z/2);
[X,Y,Z] = meshgrid((2:X_window-1).*voxel_X-r2,(2:Y_window-1).*voxel_Y-r1,(2:Z_window-1).*voxel_Z-r3);

% Initial Centers for 5 protocols
I_C = zeros([size(X),1]);


for i=1:size(X,1)
    for j = 1:size(X,2)
        for l = 1:size(X,3)
            I_C (i,j,l) = I(X(i,1,1), Y(1,j,1),Y(1,1,l));
        end
    end
end


C = [X(:),Y(:),Z(:)];
for P = 1:1
    I_C_temp = squeeze(I_C(:,:,:,P));
    C = cat(2,C,I_C_temp(:));
end

k = size(C,1);

% Update supervoxel centers
for i = 1:k
    temp_C = C(i,:);
    x = temp_C(1);
    y = temp_C(2);
    z = temp_C(3);
    
    temp_W = I_temp(x-1:x+1,y-1:y+1,z);
    a = abs(gradient(temp_W));
    [x,y] = find(a==min(a(:)), 1, 'first');
    C(i,1) = temp_C(1)+x-2;
    C(i,2) = temp_C(2)+y-2;
    
    for P = 1:1
        C(i,3+1) = I(C(i,1), C(i,2),P);
    end
    
    
end

% Assigning the lables
Label = -1*ones(size(I_temp));
Distance = inf*ones(size(I_temp));

counter = Iterations; % Iterastions
while counter
    for c_k = 1:k
        Center = [C(c_k,1),C(c_k,2),C(c_k,3)];
        for i = Center(1)-2*r1+1:Center(1)+2*r1-1
            for j = Center(2)-2*r2+1:Center(2)+2*r2-1
                for l = Center(3)-2*r3+1:Center(3)+2*r3-1
                    % Calculating the spatial distance
                    d_c_term = 0;
                    d_c_temp = 0;
                    
                    for P = 1:1
                        d_c_temp = (I(Center(1),Center(2),Center(3),P)-I(i,j,l,P))^2;
                        d_c_term = d_c_term + d_c_temp;
                    end
                    
                    d_c = sqrt(d_c_term );
                    
                    % Calculating the intensity distance
                    d_s = sqrt((term_x*(Center(1)-i))^2 + (term_y*(Center(2)-j))^2+(term_z*(Center(3)-l))^2);
                    
                    % Calculating the total distance D between C_k and i
                    D = sqrt((d_c)^2+(d_s/S)^2*Compactness^2);
                    if D < Distance(i,j,l)
                        Distance(i,j,l) = D;
                        Label(i,j,l) = c_k;
                    end
                end % for l
            end % for j
        end % for i
    end % for SV Label
    
    %% Updating stage
    counter = counter-1;
    
    % Compute new cluster centers
    if counter > 0
        for c_k = 1:k
            [x1,y1,z1] = find3((Label==c_k));
            C(c_k,1) = round(mean(x1));
            C(c_k,2) = round(mean(y1));
            C(c_k,3) = round(mean(z1));
        end
        
        % make sure that the centers are inside image dimensions
        C(C(:,1)<2*r1,1) = 2*r1;
        C(C(:,1)>size(I_temp,1)-2*r1,1) =size(I_temp,1)-2*r1;
        C(C(:,2)<2*r2,2) = 2*r2;
        C(C(:,2)>size(I_temp,2)-2*r2,2) =size(I_temp,2)-2*r2;
        C(C(:,3)<2*r3,3) = 2*r3;
        C(C(:,3)>size(I_temp,3)-2*r3,3) =size(I_temp,3)-2*r3;
    end % update
end % While
SLIC_Labels_3D = Label;
end

%% Save
% Save the supervoxel map volumes into MAT file
% Output_Name = fullfile(Output_Path,['MRI_SLIC_Labels_Size',num2str(voxel_X),...
%     'x',num2str(voxel_Y),'x',num2str(voxel_Z),'_Compactness_0',Cmpt,'_Case_',num2str(Case),'.mat']);
% save (Output_Name,'SLIC_Labels_3D');


