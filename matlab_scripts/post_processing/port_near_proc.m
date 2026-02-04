%% processing of the nearfield between two distinct nozzles
clear;clc;close all;
%% import the q_criterion version for now
indir = "/Volumes/T7 Shield/ONR_proj_temp_result/port/217M/pcprobes_near_port/"
xyz_fname = "near_port_q_only.pxyz"
data_fname = "near_port_q_only.00750000.pcd"

outdir = indir;
Nx = 541;
Ny = 229;
Nz = 200;
Npts = Nx*Ny*Nz;
%% importing data and reorganizing
xyz = importdata(indir+xyz_fname).data;
index = (xyz(:,1))+1;
x = zeros(Npts,1);
y = zeros(Npts,1);
z = zeros(Npts,1);
x(index) = xyz(:,2);
y(index) = xyz(:,3);
z(index) = xyz(:,4);
clear("xyz")
%% checks
figure; plot(x);
figure; plot(y);
figure; plot(z);
%% read in reference geometry
fv = stlread("nozzle_geometry.stl");
vertices = fv.Points;
factor = 1/0.475;
scaledVertices = vertices*factor; 
%% read in data
data_r = importdata(indir+data_fname).data; %this is only q criterion
data = zeros(Npts,1);
data(index) = data_r;
%%
X = reshape(x,Nx,Nz,Ny);
Y = reshape(y,Nx,Nz,Ny);
Z = reshape(z,Nx,Nz,Ny);
q_cri = reshape(data,Nx,Nz,Ny);
%clear("x","y","z")
%% 
figure; plot(squeeze(X(:,1,1)))
figure; plot(squeeze(Z(1,:,1)))
figure; plot(squeeze(Y(1,1,:)))
%% plot q_criterion
iso_q = 800;
figure; hold on;
p = patch(isosurface(X(:,2:end,:),Y(:,2:end,:),Z(:,2:end,:),q_cri(:,2:end,:),iso_q));
set(p, 'FaceColor','red','EdgeColor','none')
daspect([1,1,1])
view([-1,-1,1]);camlight; lighting gouraud;
%trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), scaledVertices(:, 3),'FaceColor','b','FaceAlpha',0.3)
saveas(gcf,outdir+sprintf("near_q_iso_%d.png",iso_q))
hold off;