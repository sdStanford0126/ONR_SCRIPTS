%% port injection velocity profiles - Outlet plane and port planes
%this script is based on previous scripts
%outlet_mass_flow and port mass flow
%its purpose it to 
% a: determine injection ratio
% b: plot of velocity profiles for injection ports
% c: provide outlet plane assement between baseline and injection cases
clear;clc;close all; setPlotPref(1.4,'latex',14); format long

%% section a: determine injection ratio
%% section a1: get the outlet plane data
% I/O
indir = "/Volumes/T7 Shield/ONR_proj_result/port/252M/" %on MacOS, may need to shift to Linux as needed
outdir = indir; %store at same location for now
dt        = 2.0e-4;
%% flux base outlet plane data
outFluxdir = indir + "/fluxprobes_out/"
outlet_fname = outFluxdir + "/outlet.fp"

outlet_flux_data = importdata(outlet_fname).data;

%only extract mass flux time history for now

MF_out = outlet_flux_data(:,7);

[MF_out_std, MF_out_mean] = std(MF_out);

%% extrat ports data
topFluxdir = indir + "/fluxprobes_top/";
botFluxdir = indir + "/fluxprobes_bot/";

topDatafmt = topFluxdir + "/top%d.fp"
botDatafmt = botFluxdir + "/bot%d.fp"
MF_port_mean_sum = 0;
MF_port_std_sum  = 0;
count_validation = 0
for i = 1:6
    %6 ports top and bottom
    topDataName = sprintf(topDatafmt,i)
    botDataName = sprintf(botDatafmt,i)
    
    top_flux_data = importdata(topDataName).data;
    bot_flux_data = importdata(botDataName).data;

    MF_top = top_flux_data(:,7)
    MF_bot = bot_flux_data(:,7)

    [MF_top_std, MF_top_mean] = std(MF_top);
    [MF_bot_std, MF_bot_mean] = std(MF_bot);
    
    MF_port_mean_sum = MF_port_mean_sum + MF_bot_mean + MF_top_mean;
    MF_port_std_sum  = MF_port_std_sum + MF_top_std + MF_bot_std;
    count_validation = count_validation + 2
end

MFR_mean = MF_port_mean_sum/(MF_out_mean-MF_port_mean_sum) * 100
MFR_max  = (MF_port_mean_sum+MF_port_std_sum)/(MF_out_mean - MF_out_std - MF_port_mean_sum - MF_port_std_sum) * 100
MFR_min  = (MF_port_mean_sum-MF_port_std_sum)/(MF_out_mean + MF_out_std - MF_port_mean_sum + MF_port_std_sum) * 100

%% section b: get veloicty profile for the ports


portProbedir = indir + "/pcprobes_ports_in/"

Nx = 25;
Ny = 45;
Nz = Nx; 

%% 
%automatically find the tid_start and tid_end
files = dir(portProbedir+'ports.*.pcd');
names = {files.name};

t = cellfun(@(s) ...
    str2double(regexp(s, 'ports\.(\d+)\.pcd', 'tokens', 'once')), ...
    names);
t = sort(t);
tid_start = t(1);
tid_end = t(end);
step = t(2) - t(1);

%%
%tid_start = ;
%tid_end   = ;
%step      = 600;


tid = tid_start:step:tid_end;


prefix = "ports";
suffix = ".pcd";
dataname_fmt =  prefix + ".%08d" + suffix;

pts_file = portProbedir+ prefix+".pxyz"

xyz = importdata(pts_file).data;
ind = xyz(:,1)+1;

x(ind) = xyz(:,2);
y(ind) = xyz(:,3);
z(ind) = xyz(:,4);

clear xyz;

%% ingest necessary data and organize by ports
%need all three velocities, density and pressure
Ntot  = length(ind); %number of pts for all ports
Nport = Nx*Ny;       %number of pts per port
Ntid  = length(tid); %time history length

u_t   = zeros(Ntid,Ntot);
v_t   = zeros(Ntid,Ntot);
w_t   = zeros(Ntid,Ntot);
rho_t = zeros(Ntid,Ntot);
p_t   = zeros(Ntid,Ntot);

for i = 1:length(tid)
     dataname = sprintf(dataname_fmt,tid(i));
     fprintf(dataname + "\n")
     probedata = importdata(portProbedir +dataname).data;
     u_t(i,ind) = probedata(:,1);
     v_t(i,ind) = probedata(:,2);
     w_t(i,ind) = probedata(:,3);
     p_t(i,ind) = probedata(:,4);
     rho_t(i,ind)=probedata(:,5);
end


x_shift = 0.5*3.2/2^11;
x_str = -0.28027368+x_shift;
x_end = -0.22498947-x_shift;
x_vec = linspace(x_str,x_end,Nx);

dy = 0.09212632 - 2*x_shift;
y_vec = linspace(0,dy,Ny);

%%
veln_xc_str = zeros(Nx,12);
veln_yc_str = zeros(Ny,12);

for i = 1:12
    %iterate through all ports
    %top 1-6
    %bot 7-12
    if i <7
        %for top ports
        vec1 = transpose([-0.2250-(-0.2803),0,0.4833-0.4792]);
        vec2 = transpose([0,1,0]);
    else
        %for bot ports
        vec1 = transpose([0.2250-0.2803,0,0.4833-0.4792]);
        vec2 = transpose([0,1,0]);
    end
    vec_n = cross(vec2,vec1);
    vec_n = vec_n/norm(vec_n); %3 components of normal vectors

    %get necessary indices
    ind_cur_str = (i-1)*Nport+1;
    ind_cur_end = i*Nport;
    inds = ind_cur_str:ind_cur_end;

    %get the current grid
    %x_c = x(inds);
    %y_c = y(inds);
    %z_c = z(inds);
    
    %the y centerline
    % x is center --> 13
    x_center_inds = ind_cur_str+12:Nx:ind_cur_end;
    u_yc = u_t(:,x_center_inds);
    v_yc = v_t(:,x_center_inds);
    w_yc = w_t(:,x_center_inds);

    %vel_yc = zeros(Ntid,Ny,3);

    [u_std_yc, u_mean_yc] = std(u_yc,1,1);
    [v_std_yc, v_mean_yc] = std(v_yc,1,1);
    [w_std_yc, w_mean_yc] = std(w_yc,1,1);

    vel_mean_yc = zeros(Ny,3);
    vel_mean_yc(:,1) = u_mean_yc;
    vel_mean_yc(:,2) = v_mean_yc;
    vel_mean_yc(:,3) = w_mean_yc;

    %vel_std_yc = zeros(Ny)

    veln_mean_yc = dot(vel_mean_yc,repmat(vec_n',Ny,1),2);
    veln_yc_str(:,i) = veln_mean_yc;
    %the x centerline
    % y is center --> 23
    y_center_inds = ind_cur_str+22*Nx:ind_cur_str+23*Nx-1;
    % plot diagnostics
    % figure(1)
    % scatter3(x(x_center_inds),y(x_center_inds),z(x_center_inds),'r.')
    % hold on;
    % figure(1)
    % scatter3(x(y_center_inds),y(y_center_inds),z(y_center_inds),'b.')
    % axis equal     
    u_xc = u_t(:,y_center_inds);
    v_xc = v_t(:,y_center_inds);
    w_xc = w_t(:,y_center_inds);

    [u_std_xc, u_mean_xc] = std(u_xc,1,1);
    [v_std_xc, v_mean_xc] = std(v_xc,1,1);
    [w_std_xc, w_mean_xc] = std(w_xc,1,1);

    vel_mean_xc = zeros(Nx,3);
    vel_mean_xc(:,1) = u_mean_xc;
    vel_mean_xc(:,2) = v_mean_xc;
    vel_mean_xc(:,3) = w_mean_xc;

    veln_mean_xc = dot(vel_mean_xc,repmat(vec_n',Nx,1),2);
    veln_xc_str(:,i) = veln_mean_xc;

end
%% 

[veln_xc_str_std,veln_xc_str_mean] = std(veln_xc_str,1,2);
[veln_yc_str_std,veln_yc_str_mean] = std(veln_yc_str,1,2);

veln_yc_str_mean = veln_yc_str_mean';
veln_xc_str_mean = veln_xc_str_mean';
veln_yc_str_std  = veln_yc_str_std';
veln_xc_str_std = veln_xc_str_std';
%% section b: profiles along centerline
figure;
hold on; box on;
plot(x_vec, veln_xc_str_mean,'r-')
plot(x_vec, veln_xc_str_mean+veln_xc_str_std,'r-.')
plot(x_vec, veln_xc_str_mean-veln_xc_str_std,'r-.')
xlabel("$x/h$")
ylabel("$|\vec{u}_n|/c_{\infty}$")
saveas(gcf,outdir + "port_avg_profile_xline.png")

figure;
hold on; box on;
plot(y_vec, veln_yc_str_mean,'b-')
plot(y_vec, veln_yc_str_mean+veln_yc_str_std,'b-.')
plot(y_vec, veln_yc_str_mean-veln_yc_str_std,'b-.')
xlabel("$y/h$")
ylabel("$|\vec{u}_n|/c_{\infty}$")
xlim([0,dy])
saveas(gcf,outdir + "port_avg_profile_yline.png")

save(outdir + "port_centerline_profiles.mat", ...
    "x_vec","y_vec","veln_xc_str_mean","veln_xc_str_std",...
    "veln_yc_str_mean","veln_yc_str_std")
%% section c: plotting outlet faces


prefix = "outlet";
suffix = ".pcd";
dataname_fmt = prefix + ".%08d" + suffix;

outProbedir = indir + "/pcprobes_outlet_in/"
%outProbedir = "/Volumes/T7 Shield/ONR_proj_temp_result/port/217M/pcprobes_outlet_in/"
tid_start = 750;
tid_end   = 550500;
step      = 750;
dt        = 2.0e-4;

tid = tid_start:step:tid_end;
pts_file =  outProbedir  +prefix+".pxyz"

xyz = importdata(pts_file).data;
ind = xyz(:,1)+1;

x(ind) = xyz(:,2);
y(ind) = xyz(:,3);
z(ind) = xyz(:,4);

Ny = 159;
Nz = 79;

Ntid = length(tid);
Ntot = length(ind);
delta=3.2/256;
zs = 0:delta:0.49814422;
zs = [-flip(zs(2:end)),zs];

ys = -1+delta:delta:1-delta;
u_t   = zeros(Ntid,Ntot);
v_t   = zeros(Ntid,Ntot);
w_t   = zeros(Ntid,Ntot);
rho_t = zeros(Ntid,Ntot);
p_t   = zeros(Ntid,Ntot);

for i = 1:length(tid)
     dataname = sprintf(dataname_fmt,tid(i));
     fprintf(dataname + "\n")
     probedata = importdata(outProbedir +dataname).data;
     u_t(i,ind) = probedata(:,1);
     v_t(i,ind) = probedata(:,2);
     w_t(i,ind) = probedata(:,3);
     p_t(i,ind) = probedata(:,4);
     rho_t(i,ind) = probedata(:,5);
end

%% reform the data
%plot velocity profile
port_u_mean = mean(u_t,1);
port_u_mean = reshape(port_u_mean,[Nz,Ny]);
figure; hold on; box on;
s1 = pcolor(ys,zs,port_u_mean);
s1.EdgeColor="none";
s1.FaceColor="interp";
xlabel("$y$")
ylabel("$z$")
cb = colorbar; clim([0,1.4]); axis equal;ylim([-0.5,0.5]);xlim([-1,1]);
ylabel(cb,"$\frac{u}{c_{\infty}}$",interpreter='latex',FontSize=16,Rotation=0)
saveas(gcf, outdir+"port_outlet_u_profile.png")

%plot pressure profile
port_p_mean = mean(p_t,1);
port_p_mean = reshape(port_p_mean,[Nz,Ny]);
figure; hold on; box on;
s2 = pcolor(ys,zs,port_p_mean);
colormap cool
s2.EdgeColor="none";
s2.FaceColor="interp";
xlabel("$y$")
ylabel("$z$")
cb = colorbar; clim([0,1.4]); axis equal;ylim([-0.5,0.5]);xlim([-1,1]);
ylabel(cb,"$\frac{p}{\rho_{\infty}c_{\infty}^2}$",interpreter='latex',FontSize=16,Rotation=0)
saveas(gcf, outdir+"port_outlet_p_profile.png")
%plot mean density
port_rho_mean = mean(rho_t,1);
port_rho_mean = reshape(port_rho_mean,[Nz,Ny]);
figure;hold on; box on;
s3 = pcolor(ys,zs,port_rho_mean);
s3.EdgeColor="none";
s3.FaceColor="interp";
xlabel("$y$")
ylabel("$z$")
cb = colorbar;  axis equal;ylim([-0.5,0.5]);xlim([-1,1]);
ylabel(cb,"$\frac{rho}{\rho_{\infty}}$",interpreter='latex',FontSize=16,Rotation=0)
saveas(gcf, outdir+"port_outlet_rho_profile.png")
save(outdir + "port_outlet_profile.mat","ys","zs",...
    "port_u_mean","port_p_mean","port_rho_mean")

figure;
port_v_mean = reshape(mean(v_t,1),[Nz,Ny]);
port_w_mean = reshape(mean(w_t,1),[Nz,Ny]);
quiver(ys,zs,port_v_mean,port_w_mean)

%% store port time record for turbulence data
ind_tgt = find(abs(y-0.143) < delta/2 & abs(z - min(z)-delta) < delta/2); %finding the center injection plane near the bottom layer
y_tgt = y(ind_tgt)
z_tgt = z(ind_tgt)
figure; hold on;
s2 = pcolor(ys,zs,port_u_mean);
s2.EdgeColor="none";
s2.FaceColor="interp";
plot(y_tgt,z_tgt,'r*');axis equal
saveas(gcf,"porbe_loc_port.png")

u_port_tgt = u_t(:,ind_tgt);

uf_port_tgt = u_port_tgt - mean(u_port_tgt);
L = length(uf_port_tgt);
Fs = 1/(step*dt)
f_port = Fs/L *(0:(floor(L/2)));
P2 = 1/L * fft(uf_port_tgt);
PSD1 = abs(P2(1:floor(L/2)+1)).^2;
PSD1(2:end-1) = 2*PSD1(2:end-1);
PSDport = PSD1;
figure; loglog(f_port,abs(P2(1:L/2+1)))
%%

prefix = "outlet";
suffix = ".pcd";
dataname_fmt = prefix + ".%08d" + suffix;

outProbedir = indir + "/baseline/pcprobes_outlet_in/"
tid_start = 1047150;
tid_end   = 1122000;
step      = 150;
dt        = 1e-3;

tid = tid_start:step:tid_end;
pts_file =  outProbedir  +prefix+".pxyz"

xyz = importdata(pts_file).data;
ind = xyz(:,1)+1;

x(ind) = xyz(:,2);
y(ind) = xyz(:,3);
z(ind) = xyz(:,4);

Ny = 159;
Nz = 79;

Ntid = length(tid);
Ntot = length(ind);

zs = 0:delta:0.49814422;
zs = [-flip(zs(2:end)),zs];

ys = -1+delta:delta:1-delta;
u_t   = zeros(Ntid,Ntot);
%v_t   = zeros(Ntid,Ntot);
%w_t   = zeros(Ntid,Ntot);
rho_t = zeros(Ntid,Ntot);
p_t   = zeros(Ntid,Ntot);

for i = 1:length(tid)
     dataname = sprintf(dataname_fmt,tid(i));
     fprintf(dataname + "\n")
     probedata = importdata(outProbedir +dataname).data;
     u_t(i,ind) = probedata(:,1);
     %v_t(i,ind) = probedata(:,2);
     %w_t(i,ind) = probedata(:,3);
     p_t(i,ind) = probedata(:,4);
     rho_t(i,ind) = probedata(:,5);
end

%% reform the data
baseline_u_mean = mean(u_t,1);
baseline_u_mean = reshape(baseline_u_mean,[Nz,Ny]);
figure; hold on; box on;
s1 = pcolor(ys,zs,baseline_u_mean);
s1.EdgeColor="none";
s1.FaceColor="interp";
xlabel("$y$")
ylabel("$z$")
cb = colorbar; clim([0,1.4]); axis equal;ylim([-0.5,0.5]);xlim([-1,1]);
ylabel(cb,"$\frac{u}{c_{\infty}}$",interpreter='latex',FontSize=16,Rotation=0)
saveas(gcf, outdir+"baseline_outlet_u_profile.png")




%plot pressure profile
baseline_p_mean = mean(p_t,1);
baseline_p_mean = reshape(baseline_p_mean,[Nz,Ny]);
figure; hold on; box on;
s2 = pcolor(ys,zs,baseline_p_mean);
colormap cool
s2.EdgeColor="none";
s2.FaceColor="interp";
xlabel("$y$")
ylabel("$z$")
cb = colorbar; clim([0,1.4]); axis equal;ylim([-0.5,0.5]);xlim([-1,1]);
ylabel(cb,"$\frac{p}{\rho_{\infty}c_{\infty}^2}$",interpreter='latex',FontSize=16,Rotation=0)
saveas(gcf, outdir+"baseline_outlet_p_profile.png")


%% store port time record for turbulence data
ind_tgt = find(abs(y-0.143) < delta/2 & abs(z - min(z)-delta) < delta/2); %finding the center injection plane near the bottom layer
y_tgt = y(ind_tgt)
z_tgt = z(ind_tgt)
figure; hold on;
s2 = pcolor(ys,zs,baseline_u_mean);
s2.EdgeColor="none";
s2.FaceColor="interp";
plot(y_tgt,z_tgt,'r*')   
u_base_tgt = u_t(:,ind_tgt); axis equal
saveas(gcf,"porbe_loc_base.png")
uf_base_tgt = u_base_tgt - mean(u_base_tgt);
L = length(uf_base_tgt);
Fs = 1/(step*dt)
f_base = Fs/L *(0:(floor(L/2)));
P2 = 1/L * fft(uf_base_tgt);
PSD1 = abs(P2(1:floor(L/2)+1)).^2;
PSD1(2:end-1) = 2*PSD1(2:end-1);
PSDbase = PSD1;
figure; loglog(f_base,abs(P2(1:L/2+1)))

%% plot single point comparisons
figure; 
loglog(f_port,PSDport,'r-',DisplayName="port, 217M")
hold on
loglog(f_base,PSDbase,'b:',DisplayName="Baseline, 144M")
xlabel("$f/(c_\infty/h)$")
ylabel("$PSD\: u'$")
title("$outlet\: u'\: PSD comp.$",FontWeight="bold")
legend(Location="southwest")
saveas(gcf,"up_comp_PSD.png")


%% Build thrust profile
%get baseline density profile (also plot just for reference)
baseline_rho_mean = mean(rho_t,1);
baseline_rho_mean = reshape(baseline_rho_mean,[Nz,Ny]);
figure; hold on; box on;
s3 = pcolor(ys,zs,baseline_rho_mean);
s3.EdgeColor="none";
s3.FaceColor="interp";
cb = colorbar; axis equal;ylim([-0.5,0.5]);xlim([-1,1]);
ylabel(cb,"$\frac{rho}{\rho_{\infty}}$",interpreter='latex',FontSize=16,Rotation=0)
saveas(gcf, outdir+"baseline_outlet_rho_profile.png")
save(outdir + "baseline_outlet_profile.mat","ys","zs",...
    "baseline_u_mean","baseline_p_mean","baseline_rho_mean")
%ambient pressure
gamma = 1.4;
p_inf = 1/gamma;


%compute massflow
Ly = max(ys) - min(ys)
Lz = max(zs) - min(zs)
A = Ly*Lz %area of the acutal sample nozzle due to size reduction to avoid 
% missed sampled points

MF_baseline = trapz(zs,trapz(ys, baseline_rho_mean .* baseline_u_mean,2))
MF_port     = trapz(zs,trapz(ys, port_rho_mean .* port_u_mean,2))








%% Thurst compute
baseline_mom_flux = baseline_u_mean.^2.* baseline_rho_mean + ...
    baseline_p_mean - p_inf;
port_mom_flux     = port_u_mean.^2.*port_rho_mean + port_p_mean - p_inf;

F_thr_baseline = trapz(zs,trapz(ys,baseline_mom_flux,2))
F_thr_port     = trapz(zs,trapz(ys,port_mom_flux,2))

%thrust to mass flow ratio
T2MFR_baseline = F_thr_baseline/MF_baseline
T2MFR_port = F_thr_port/MF_port

F_thr_ideal = MF_port*T2MFR_baseline

F_thr_lost = ( F_thr_ideal- F_thr_port)/F_thr_ideal*100