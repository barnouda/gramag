function z = demalt(topofilename,x,y)

[X,Y,Z] = readTopoFile(topofilename);
z = interp2(X,Y,Z,x,y);

