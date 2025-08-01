%%% Vérif example2d.c (cubature)

clear all
close all

xmin = -2;
xmax = 2;
ymin = -2;
ymax = 2;

Q = dblquad(@f,xmin,xmax,ymin,ymax,1.0e-4);


