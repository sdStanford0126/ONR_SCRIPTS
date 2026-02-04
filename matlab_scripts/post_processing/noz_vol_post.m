%processing for the nozzle interior probes
clear;clc;close all
%% read the data
indir = "/Volumes/T7 Shield/ONR_proj_temp_result/port/217M/pcprobes_noz_vol/"
xyz_fname = "noz_vol.pxyz"
data_fname = "noz_vol.00750000.pcd"

outdir = indir;

%pxyz data format are index, x,y,z
xyz = importdata(indir+xyz_fname).data;
Nx_noz = 81;
Ny_noz = 159;
Nz_noz = 79;
%%
Npt = size(xyz,1);
index = xyz(:,1)+1;
x = zeros(Npt,1);
y = zeros(Npt,1);
z = zeros(Npt,1);

x(index)=xyz(:,2);
y(index)=xyz(:,3);
z(index)=xyz(:,4);
%% 
fv = stlread("nozzle_geometry.stl");
vertices = fv.Points;
factor = 1/0.475;
scaledVertices = vertices*factor; 
%% quick check on reorgianzation
% figure; hold on;
% plot(x,DisplayName="x")
% legend(); hold off;
% figure; hold on;
% plot(y,DisplayName="y")
% legend(); hold off;
% figure; hold on;
% plot(z,DisplayName="z")
% legend(); hold off;
%close all;
%% load data file
data_r = importdata(indir+data_fname).data;
% data labels are as follow
%1 avg(p), 2 avg(rho), 3 avg(T), 4 avg(u), 5 avg(v), 6 avg(w),
% 7 rms(u), 8 rms(v), 9 rms(w), 10 u, 11 v, 12 w,
%13 lambda2, 14 mag(avg(grad(p)), 15 mag(avg(grad(rho)), 
% 16 mag(grad(p)), 17 mag(grad(rho)),18 p, 19 q_criterion, 20 rho, 21 T

%avg flow datas
p_avg   = zeros(Npt,1);
rho_avg = zeros(Npt,1);
u_avg   = zeros(Npt,1);
v_avg   = zeros(Npt,1);
w_avg   = zeros(Npt,1);

p_avg(index) = data_r(:,1);
rho_avg(index) = data_r(:,2);
u_avg(index) = data_r(:,4);
v_avg(index) = data_r(:,5);
w_avg(index) = data_r(:,6);

%q_criterion visualization
q_cri = zeros(Npt,1);
q_cri(index) = data_r(:,19);

%density gradient
grad_rho_mag = zeros(Npt,1);
grad_rho_mag(index) = data_r(:,15);

%% reshape the arrays
X = reshape(x,Nz_noz,Ny_noz,Nx_noz);
Y = reshape(y,Nz_noz,Ny_noz,Nx_noz);
Z = reshape(z,Nz_noz,Ny_noz,Nx_noz);
% figure; hold on;
% plot(squeeze(X(1,1,:)),DisplayName="x")
% plot(squeeze(Y(1,:,1)),DisplayName="y")
% plot(squeeze(Z(:,1,1)),DisplayName="z")
% legend; hold off;



%% quick q_criterion visualization
iso_q = 10;
figure; hold on;
p = patch(isosurface(X,Y,Z,reshape(q_cri,Nz_noz,Ny_noz,Nx_noz),iso_q));
set(p, 'FaceColor','red','EdgeColor','none')
daspect([1,1,1])
view([-1,-1,1]);camlight; lighting gouraud;
%trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), scaledVertices(:, 3),'FaceColor','b','FaceAlpha',0.3)
saveas(gcf,outdir+"noz_q_iso_10.png")
hold off;
%% u_avg quick visualization
mach = sqrt(u_avg.^2 + v_avg.^2 + w_avg.^2);
iso_mach = 0.99;
figure;hold on;
p = patch(isosurface(X,Y,Z,reshape(mach,Nz_noz,Ny_noz,Nx_noz),iso_mach));
set(p, 'FaceColor','red','EdgeColor','none')
daspect([1,1,1])
view([-1,-1,1]);camlight; lighting gouraud;
hold off;
%% grad_rho_mag quick visualization
iso_grad_rho_mag = 5;
Grad_rho = reshape(grad_rho_mag,Nz_noz,Ny_noz,Nx_noz);

figure;hold on;
p = patch(isosurface(X(2:end-1,2:end-1,2:end-1),Y(2:end-1,2:end-1,2:end-1)...
    ,Z(2:end-1,2:end-1,2:end-1),Grad_rho(2:end-1,2:end-1,2:end-1),iso_grad_rho_mag));
set(p, 'FaceColor','red','EdgeColor','none')
daspect([1,1,1])
view([-1,-1,0.5]);camlight; lighting gouraud;
hold off;
saveas(gcf,outdir + "noz_shock_view1.png")
view([1,-1,0.5])
saveas(gcf,outdir+"noz_shock_view2.png")
saveas(gcf,outdir+"noz_shock.fig")