function [XI,YI,ZI] = readTopoFile(file)
fp = fopen(file,'r');
latmin = fread(fp,1,'float32');
latmax = fread(fp,1,'float32');
lonmin = fread(fp,1,'float32');
lonmax = fread(fp,1,'float32');
DXY = fread(fp,1,'float32');
dx = DXY;
dy = DXY;
xi = lonmin:dx:lonmax;
yi = latmin:dy:latmax;
[XI,YI] = meshgrid(xi,yi);
n = fread(fp,1,'int32');
m = fread(fp,1,'int32');
ZI = fread(fp,[n,m],'float32');
ZI = ZI';
fclose(fp);