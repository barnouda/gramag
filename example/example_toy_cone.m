%
% Synthetic model of a toy volcanic cone
% with checkerboard density distribution
%

clear all
close all

show_figures = 1;
write_outputs = 1; % write input files for modelling

%% Topography: build and write topography binary file

xmin_topo = -2000;
xmax_topo = 2000;
ymin_topo = -2000;
ymax_topo = 2000;

dxyz_topo = 10;

x_topo = xmin_topo:dxyz_topo:xmax_topo;
y_topo = ymin_topo:dxyz_topo:ymax_topo;
[Y_topo,X_topo] = meshgrid(y_topo,x_topo);

Z_topo = zeros(size(X_topo));

cone_radius = 400;
cone_summit = 400;

imin = find(x_topo==-cone_radius);
imax = find(x_topo==cone_radius);
jmin = find(y_topo==-cone_radius);
jmax = find(y_topo==cone_radius);

for i = imin:imax
    for j = jmin:jmax
        s = sqrt(x_topo(i)^2+y_topo(j)^2);
        if s <= cone_radius
            Z_topo(j,i) = cone_summit-s;
        end
    end
end

if write_outputs == 1
    writeTopoFile('toy_cone_topo.bin',ymin_topo,ymax_topo,xmin_topo,xmax_topo,dxyz_topo,Z_topo);
end

%% Grid of nodes

xmin_grid = -1500;
xmax_grid = 1500;
ymin_grid = -1500;
ymax_grid = 1500;
zmin_grid = -400;
zmax_grid = 800;

dxyz_grid = 100;
x_grid = xmin_grid:dxyz_grid:xmax_grid;
y_grid = ymin_grid:dxyz_grid:ymax_grid;
z_grid = zmin_grid:dxyz_grid:zmax_grid;

[Y_grid,X_grid,Z_grid] = meshgrid(y_grid,x_grid,z_grid);
xx_grid = reshape(X_grid,numel(X_grid),1);
yy_grid = reshape(Y_grid,numel(Y_grid),1);
zz_grid = reshape(Z_grid,numel(Z_grid),1);


%% Paramters: density and magnetisation

i = 1:length(x_grid);
j = 1:length(y_grid);
k = 1:length(z_grid);
[I,J,K] = meshgrid(i,j,k);
ii = reshape(I,numel(I),1);
jj = reshape(J,numel(J),1);
kk = reshape(K,numel(K),1);

% Checkerboard of density
dd = 500; % size of checkerboard element
dr = 200; % density difference
rho = dr-2*dr*mod(floor(xx_grid/dd)+floor(yy_grid/dd),2);

% Bands of positive and negative magnetisations

dd = 1000; % size of bands
dm = 3; % magnetisation difference
mag = dm-2*dm*mod(floor(yy_grid/dd),2);

% Set parameters to zero above topography
zz_topo = interp2(Y_topo,X_topo,Z_topo,yy_grid,xx_grid);
for k = 1:length(zz_grid)
    if -zz_grid(k)>zz_topo(k)
        rho(k) = 0;
        mag(k) = 0;
    end
end


if write_outputs == 1

    % Density grid
    
    name = 'toy_cone_density_grid';

    file = fopen([name '.bin'],'w');
    fwrite(file,rho,'single');
    fclose(file);

    file = fopen([name '.par'],'w');
    fprintf(file,'%d\n',1);
    fprintf(file,'%.6f %.6f %.6f\n',xmin_grid,ymin_grid,zmin_grid);
    fprintf(file,'%d %d %d\n',length(x_grid),length(y_grid),length(z_grid));
    fprintf(file,'%.6f %.6f %.6f\n',dxyz_grid,dxyz_grid,dxyz_grid);
    fclose(file);
    
    % Magnetisation grid

    name = 'toy_cone_magnetisation_grid';

    file = fopen([name '.bin'],'w');
    fwrite(file,rho,'single');
    fclose(file);

    file = fopen([name '.par'],'w');
    fprintf(file,'%d\n',1);
    fprintf(file,'%.6f %.6f %.6f\n',xmin_grid,ymin_grid,zmin_grid);
    fprintf(file,'%d %d %d\n',length(x_grid),length(y_grid),length(z_grid));
    fprintf(file,'%.6f %.6f %.6f\n',dxyz_grid,dxyz_grid,dxyz_grid);
    fclose(file);
    
end


%% Données

xg = [];
yg = [];
ng = 0;

for xxx = -1100:200:1100
    for yyy = -1100:200:1100
        ng = ng+1;
        xg(ng) = xxx;
        yg(ng) = yyy;
    end
end

zg = demalt('toy_cone_topo.bin',xg,yg);

if write_outputs == 1
    file = fopen('toy_cone_geometry_gravi.txt','w');
    for k = 1:ng
        fprintf(file,'%.3f %.3f %.3f\n',xg(k),yg(k),zg(k));
    end
    fclose(file);
end

% Muographiques

nm = 2;

dm = 1200;
xm = [-dm 0];
ym = [0 -dm];
zm = [0 0];

daz = 5;
del = 5;

azim1 = 60:daz:120;
elev1 = 5:del:20;

azim2 = -30:daz:30;
elev2 = 5:del:20;

s = 0:500:2000;

if write_outputs == 1
%     file = fopen('toy_cone_geometry_muons.txt','w');
%     fprintf(file,'%d\n',nm);
%     fprintf(file,'%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n', ...
%         xm(1),ym(1),zm(1), ...
%         min(azim1),max(azim1),daz, ...
%         min(elev1),max(elev1),del);
%     fprintf(file,'%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n', ...
%         xm(2),ym(2),zm(2), ...
%         min(azim2),max(azim2),daz, ...
%         min(elev2),max(elev2),del);
%     fclose(file);

    file = fopen('toy_cone_geometry_muog.txt','w');
    for iel = 1:length(elev1)
        el = elev1(iel);
        for iaz = 1:length(azim1)
            az = azim1(iaz);
            fprintf(file,'%.3f %.3f %.3f %.3f %.3f\n', ...
                xm(1),ym(1),zm(1),az,el);
        end
    end
    for iel = 1:length(elev2)
        el = elev2(iel);
        for iaz = 1:length(azim2)
            az = azim2(iaz);
            fprintf(file,'%.3f %.3f %.3f %.3f %.3f\n', ...
                xm(2),ym(2),zm(2),az,el);
        end
    end
    fclose(file);
    
%     file = fopen('toy_cone_data_muons_inv.txt','w');
%     for iel = 1:length(elev1)
%         for iaz = 1:length(azim1)
%             fprintf(file,'%.3f %.3f %.3f %.3f %.3f %.3f %.3f\n',);
%         end
%     end
%     fclose(file);
end

% Magnetic data

xd = [];
yd = [];
nd = 0;

for xxx = -1100:200:1100
    for yyy = -1100:200:1100
        nd = nd+1;
        xd(nd) = xxx;
        yd(nd) = yyy;
    end
end

zd = 600*ones(size(xd));

if write_outputs == 1
    
    file = fopen('toy_cone_geometry_mag.txt','w');
    for k = 1:nd
        fprintf(file,'%.3f %.3f %.3f\n',xd(k),yd(k),zd(k));
    end
    fclose(file);
    
    % File of magnetisation and geomagnetic field directions
    file = fopen('dir.par','w');
    fprintf(file, '34780 0 19758\n');
    fprintf(file, '27363 -5590 28398\n');
    fclose(file);
    
end

%% Input parameters files for modelling code

if write_outputs == 1
   
    % Input file for gravimetric modelling
    file = fopen('input_gravi_modelling.in','w');
    fprintf(file, 'toy_cone_gravi_modelling\n');
    fprintf(file, 'toy_cone_topo.bin\n');
    fprintf(file, 'toy_cone_density_grid.par\n');
    fprintf(file, 'toy_cone_density_grid.bin\n');
    fprintf(file, 'toy_cone_geometry_gravi.txt\n');
    fprintf(file, '-1800 1800\n');
    fprintf(file, '-1800 1800\n');
    fprintf(file, '-400 800\n');
    fprintf(file, '0\n');
    fprintf(file, '1 800\n');
    fclose(file);
    
    % Input file for muographic modelling (panda code)
    file = fopen('input_muog_modelling.in','w');
    fprintf(file, 'toy_cone_gravi_modelling\n');
    fprintf(file, 'toy_cone_topo.bin\n');
    fprintf(file, 'toy_cone_density_grid.par\n');
    fprintf(file, 'toy_cone_density_grid.bin\n');
    fprintf(file, 'toy_cone_geometry_muog.txt\n');
    fclose(file);
    
    % Input file for magnetic modelling
    file = fopen('input_mag_modelling.in','w');
    fprintf(file, 'toy_cone_mag_modelling\n');
    fprintf(file, 'toy_cone_topo.bin\n');
    fprintf(file, 'toy_cone_magnetisation_grid.par\n');
    fprintf(file, 'toy_cone_magnetisation_grid.bin\n');
    fprintf(file, 'toy_cone_geometry_mag.txt\n');
    fprintf(file, '-1800 1800\n');
    fprintf(file, '-1800 1800\n');
    fprintf(file, '-400 800\n');
    fprintf(file, '0\n');
    fprintf(file, '1 800\n');
    fclose(file);
    
end

%% Figures

if show_figures == 1
    
    % Density model with gravimetric and muographic data
    figure
    set(gcf,'color','w');
    % Topography
    surf(x_topo,y_topo,-Z_topo,'edgecolor','none')
    colorbar
    colormap(gray)
    caxis([-2000 1000])
    hold on
    % Grid
    % plot3(xx_grid,yy_grid,zz_grid,'ko');
    % Densities
    ind = find(rho>0);
    plot3(xx_grid(ind),yy_grid(ind),zz_grid(ind),'mo','markersize',4,'markerfacecolor','m');
    ind = find(rho<0);
    plot3(xx_grid(ind),yy_grid(ind),zz_grid(ind),'go','markersize',4,'markerfacecolor','g');
    colorbar
    % Gravimetric data
    plot3(xg,yg,-zg,'bo','markersize',4,'markerfacecolor','b');
    % Muographic data
    plot3(xm(1),ym(1),zm(1),'rs','markersize',16,'markerfacecolor','r')
    plot3(xm(2),ym(2),zm(2),'ys','markersize',16,'markerfacecolor','y')
    for iel = 1:length(elev1)
        el = elev1(iel);
        for iaz = 1:length(azim1)
            az = azim1(iaz);
            xs = xm(1)+s*sind(az)*cosd(el);
            ys = ym(1)+s*cosd(az)*cosd(el);
            zs = zm(1)+s*sind(el);
            plot3(xs,ys,-zs,'r-');
        end
    end
    for iel = 1:length(elev2)
        el = elev2(iel);
        for iaz = 1:length(azim2)
            az = azim2(iaz);
            xs = xm(2)+s*sind(az)*cosd(el);
            ys = ym(2)+s*cosd(az)*cosd(el);
            zs = zm(2)+s*sind(el);
            plot3(xs,ys,-zs,'y-');
        end
    end
    % Figure parameters
    axis equal
    set(gca,'zdir','reverse')
    title('Density nodes, gravimetric and muographic data')
    axis equal
    hold on
    
    % Magnetisation model with magnetic data
    figure
    set(gcf,'color','w');
    % Topography
    surf(x_topo,y_topo,-Z_topo,'edgecolor','none')
    colorbar
    colormap(gray)
    caxis([-2000 1000])
    hold on
    % Grid
    plot3(xx_grid,yy_grid,zz_grid,'ko');
    % Magnetisations
    ind = find(mag>0);
    plot3(xx_grid(ind),yy_grid(ind),zz_grid(ind),'mo','markersize',4,'markerfacecolor','m');
    ind = find(mag<0);
    plot3(xx_grid(ind),yy_grid(ind),zz_grid(ind),'go','markersize',4,'markerfacecolor','c');
    colorbar
    % Magnetic data
    plot3(xd,yd,-zd,'bo','markersize',4,'markerfacecolor','b');
    % Figure parameters
    axis equal
    set(gca,'zdir','reverse')
    title('Magnetisation nodes and magnetic data')
    axis equal

end

