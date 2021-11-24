function G=readrawPOPImeta(filename)
%readraw - read RAW format grey scale image file from Disk
%	Usuage  : G=readraw(filename)

	disp(['	Retrieving Image ' filename ' ...']);
	fid=fopen(filename,'rb');
	if (fid==-1)
	  	error('can not open imput image filem press CTRL-C to exit \n');
	  	pause
	end
	pixel=fread(fid,inf, 'short');
	fclose(fid);
    X=482;
    Y=360;
    Z=141;
%    [Y X]=size(pixel);
    G=reshape(pixel,[X,Y,Z]);
%    Size=(Y*X);
%    N=sqrt(Size);
%    G=zeros(N,N);
%    G(1:Size)=pixel(1:Size);
%    G=G';
   end