function writeTopoFile(file,latmin,latmax,lonmin,lonmax,DXY,ZI)
fp = fopen(file,'w');
fwrite(fp,latmin,'float32');
fwrite(fp,latmax,'float32');
fwrite(fp,lonmin,'float32');
fwrite(fp,lonmax,'float32');
fwrite(fp,DXY,'float32');
[m,n] = size(ZI);
fwrite(fp,n,'int32');
fwrite(fp,m,'int32');
fwrite(fp,ZI','float32');
fclose(fp);