
function [info,vf]=readvf(filename)
%readraw - read RAW format grey scale image file from Disk
%	Usuage  : G=readraw(filename)

	disp(['	Retrieving Image ' filename ' ...']);
	fid=fopen(filename,'rb');
    if (fid==-1)
	  	error('can not open imput image filem press CTRL-C to exit \n');
	  	pause
    end
    head=fscanf(fid,'%20s',1);
    info.head=head;
    magic_number=fscanf(fid,'%10s',1);
    info.magic_number=magic_number;
    x=fscanf(fid,'%10s',1);
    y=fscanf(fid,'%10s',1);
    z=fscanf(fid,'%10s',1);
    
    size_x=str2double(x);
    size_y=str2double(y);
    size_z=str2double(z);
    info.size=[size_x,size_y,size_z];
    
    grid_x=fscanf(fid,'%10s',1);
    grid_y=fscanf(fid,'%10s',1);
    grid_z=fscanf(fid,'%10s',1);
    info.grid=[grid_x,grid_y,grid_z];
	pixel=fread(fid,inf, 'float');
	fclose(fid);
    vf=reshape(pixel,[size_x,size_y,size_z,3]);

 end
