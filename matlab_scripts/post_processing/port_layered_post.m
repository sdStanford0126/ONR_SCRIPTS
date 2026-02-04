%% Post proc of multilayered ports (definition in /Simulation Utilities/pcprobes_volumetric.m)
clear;clc; close all;
%% Probe configuration and I/O
indir = "/Volumes/T7 Shield/ONR_proj_temp_result/port/217M/pcprobes_layered_port/";
outdir = "./";
port_top_pxyz = "top_ports.pxyz";
port_bot_pxyz = "bot_ports.pxyz";
port_top_fname = "top_ports.00750000.pcd";
port_bot_fname = "bot_ports.00750000.pcd";

Nx = 43;
Ny = 73;
Nz = 43;
numPort = 6; %number of ports on either top/bottom
Nlayers = 4; %number of layers per port
Npts = Nx*Ny*numPort*Nlayers;
%% load and reorganize data
xyz_top_r = importdata(indir + port_top_pxyz).data;
index = xyz_top_r(:,1)+1; %shift index due to matlab indexing

x = zeros(Npts,1);
y = zeros(Npts,1);
z = zeros(Npts,1);

x(index) = xyz_top_r(:,2);
y(index) = xyz_top_r(:,3);
z(index) = xyz_top_r(:,4);

%% check organization
%figure; plot(x);
%figure; plot(y);
%figure; plot(z);

%only really need to know y locations and also normalized x distance
%get the first layer of the first port and let's do this
x_port1 = x(1:Nx*Ny);
y_port1 = y(1:Nx*Ny);
z_port1 = z(1:Nx*Ny);

%%
figure; hold on;
plot(x_port1,DisplayName="x")
plot(y_port1,DisplayName="y")
plot(z_port1,DisplayName="z")
legend()

%% getting x,y,z vectors;
x_vec_top=x_port1(1:Nx);
y_vec_top=y_port1(1:Nx:end);
z_vec_top=z_port1(1:Nx);

dy_int = y_vec_top - y_port1(1);
dx_int = sqrt((x_vec_top - x_port1(1)).^2 + (z_vec_top-z_port1(1)).^2);
figure; plot(dy_int);
figure; plot(dx_int);

%plot check
%figure; hold on;
%plot(x_vec_top,DisplayName="x")
%plot(y_vec_top,DisplayName="y")
%plot(z_vec_top,DisplayName="z")
%legend()

%% for top ports
vec1 = transpose([-0.2250-(-0.2803),0,0.4833-0.4792]);
vec2 = transpose([0,1,0]);
vec_n = cross(vec2,vec1);
vec_n_top = vec_n/norm(vec_n); %3 components of normal vectors

%% read in data for top ports
data_top_r = importdata(indir+port_top_fname).data;
rho_avg = zeros(Npts,1);
vel_avg = zeros(Npts,3);
rho_avg(index) = data_top_r(:,2);
vel_avg(index,:) = data_top_r(:,4:6);

%% compute mass flux for each port and each layer
massflux_top = zeros(numPort,Nlayers);
for i = 1:numPort
    for j = 1:Nlayers
        ind_str = (i-1)*Nx*Ny*Nlayers + (j-1)*Nx*Ny+1;
        ind_end = ind_str+Nx*Ny-1;
        vel_layer = vel_avg(ind_str:ind_end,:);
        rho_layer = rho_avg(ind_str:ind_end);
        veln_layer = zeros(Nx,Ny);
        rhon_layer = zeros(Nx,Ny);
        for ix = 1:Nx
            for iy = 1:Ny
                ind = (iy-1)*Nx + ix;
                veln_layer(ix,iy) = dot(vel_layer(ind,:),vec_n_top);
                rhon_layer(ix,iy) = rho_layer(ind);
            end
        end
        mass_flux_local = trapz(dx_int,trapz(dy_int,rhon_layer.*veln_layer,2));
        massflux_top(i,j) = mass_flux_local;
    end
end

%% for bot ports
vec1 = transpose([0.2250-0.2803,0,0.4833-0.4792]);
vec_n = cross(vec2,vec1);
vec_n_bot = vec_n/norm(vec_n); %3 components of normal vectors
%% read in data for bot ports
data_bot_r = importdata(indir+port_bot_fname).data;
rho_avg = zeros(Npts,1);
vel_avg = zeros(Npts,3);
rho_avg(index) = data_bot_r(:,2);
vel_avg(index,:) = data_bot_r(:,4:6);

massflux_bot = zeros(numPort,Nlayers);
for i = 1:numPort
    for j = 1:Nlayers
        ind_str = (i-1)*Nx*Ny*Nlayers + (j-1)*Nx*Ny+1;
        ind_end = ind_str+Nx*Ny-1;
        vel_layer = vel_avg(ind_str:ind_end,:);
        rho_layer = rho_avg(ind_str:ind_end);
        veln_layer = zeros(Nx,Ny);
        rhon_layer = zeros(Nx,Ny);
        for ix = 1:Nx
            for iy = 1:Ny
                ind = (iy-1)*Nx + ix;
                veln_layer(ix,iy) = dot(vel_layer(ind,:),vec_n_bot);
                rhon_layer(ix,iy) = rho_layer(ind);
            end
        end
        mass_flux_local = trapz(dx_int,trapz(dy_int,rhon_layer.*veln_layer,2));
        massflux_bot(i,j) = mass_flux_local;
    end
end
%%  assess mass flux variation
massflux_top_max = max(max(massflux_top));
massflux_top_min = min(min(massflux_top));
massflux_bot_max = max(max(massflux_bot));
massflux_bot_min = min(min(massflux_bot));

massflux_min = min(massflux_bot_min,massflux_top_min)
massflux_max = max(massflux_top_max,massflux_bot_max)

massflux_mean = mean(mean(massflux_top))/2 + mean(mean(massflux_bot))/2
%% 
indir = "/Volumes/T7 Shield/ONR_proj_temp_result/MFR_vel_prof/" %on MacOS, may need to shift to Linux as needed
%outdir = indir; %store at same location for now

%% flux base outlet plane data
outFluxdir = indir + "/port/fluxprobes_out/"
outlet_fname = outFluxdir + "/outlet.fp"

outlet_flux_data = importdata(outlet_fname).data;

%only extract mass flux time history for now

MF_out = outlet_flux_data(:,7);

[MF_out_std, MF_out_mean] = std(MF_out);

%% MFR new assessment
MFR_mean = 12*massflux_mean/(MF_out_mean-massflux_mean)*100
MFR_max = 12*massflux_max/(MF_out_mean - MF_out_std - massflux_max)*100
MFR_min = 12*massflux_min/(MF_out_mean + MF_out_std - massflux_min)*100