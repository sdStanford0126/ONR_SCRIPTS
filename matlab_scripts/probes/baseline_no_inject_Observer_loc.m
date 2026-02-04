%% script use to setup init. observer location
%and setup parameters for FW-H solver
%last update: Nov. 4, 2024
%By Steven Dai (sd0126@stanford.edu)
format default
%clear;clc;close all
%% geometry consideration
% length nondim by h
h = 0.475 * 0.0254;%m %dimension of nozzle in physical domain
L = 2*h;
De = 1.3*(L*h)^(0.625)/(L+h)^(0.25) %this is corrected
%%
% i am going to do the following
% judging from what Kaurab has shared, I am assuming ideally isentropic
% expansion to esimate the Mach number of the outflow (not accurate)
% also I am assuming that the speed the sound in the jet flow is the same
% as that of the coflow (due to no additional heating than
% frictional/dissipative heating from turbulence.

NPR = 3; %total pressure ratio of nozzle inlet
gamma = 1.4; % specific heat ratio
Ma = jetMach_isen_ideal(NPR,gamma) %this is not computing Mj, but Ma


c_inf = 346.1 %m/s from reference value at 298K (wiki for now)
Uj = Ma*c_inf

%% 
dt = 0.001;
step = 150;
dsamp = step*dt;

atu = h/340;

f_physical = 1/(dsamp*atu)
%% experiment acoustic result configurations
bw = 50; %Hz, bandwith from experiment
%switch to nondim bandwidth is
St_bw = bw * De/Uj

%% geometry setup in FW-H
% _s denotes simulation
h_s = 1;
L_s = h_s *2;

De_s = 1.3*(L_s*h_s)^(0.625)/(L_s+h_s)^(0.25)
R_s = 100 * De_s
thetas = (40:150) * pi/180 ; %radian, converted from degrees 
thetas = [thetas,-thetas]
%%
bw = 50 %50
St_bw = bw * De/Uj
bw_s = St_bw * Ma/De_s
dt = 1e-3;
facq_s = 1/(50*dt)
blk_s = facq_s / bw_s

%% Setup command format
name_fmt = "OBSERVER NAME=%s_Mic%03d GEOM=POINT %.4f %.4f %.4f OPTIONS=P_HAT #theta=%d deg"
fprintf("\n")
test_name = sprintf(name_fmt,'Major',10,-123.30,0,1,30)

%based on shared slides from Kaurab and simulation setup (i.e. minor axis
%is y, major axis is z that): Major side is on y = 0 plane with varying x z
% ,and minor size is on z=0 plane with varying x y
%% setup for angle from -x axis in the z axis plane (plane normal to z axis) (varying in y),
% longer side of nozzle is parallel to this plane
% equivalent to along major in Jinah's input
% generate the format directly

z_plane_locs = zeros(length(thetas),4);
for i = 1:length(thetas)
    theta = thetas(i);
    %z-plane, x
    x = -1 * R_s * cos(theta); %-1 since theta is measure w.r.t -x axis
    z_plane_locs(i,1) = x;
    %z-plane, y
    y = R_s * sin(theta); % should not matter much due to planar symmetry
    z_plane_locs(i,2) = y;
    %z-plane, z
    z = 0; %in plane
    z_plane_locs(i,3) = z;
    %angle
    z_plane_locs(i,4) = theta*180/pi;
end
% write out the z_locs_in_plane
%writematrix(round(z_plane_locs,4),'z_plane_locs.txt','Delimiter',' ')

%% setup for angle from -x in the y axis plane (plane normal to y axis) varying in z
% shorter sides of the nozzle are parallel to this plane
% equivalent to along minor axis in Jinha's input
y_plane_locs = zeros(length(thetas),4);
for i = 1:length(thetas)
    theta = thetas(i);
    %y-plane, x
    x = -1 * R_s * cos(theta); %-1 since theta is measure w.r.t -x axis
    y_plane_locs(i,1) = x;
    %y-plane, y
    y = 0;
    y_plane_locs(i,2) = y;
    %y-plane, z
    z = R_s * sin(theta); % should not matter much due to planar symmetry
    y_plane_locs(i,3) = z;
    %angle
    y_plane_locs(i,4) = theta*180/pi;
end
% write out the z_locs_in_plane
%writematrix(round(y_plane_locs,4),'y_plane_locs.txt','Delimiter',' ')

figure; 
scatter3(z_plane_locs(:,1),z_plane_locs(:,2), z_plane_locs(:,3),'filled','bo'), hold on
scatter3(y_plane_locs(:,1),y_plane_locs(:,2),y_plane_locs(:,3),'filled','ro')
xlabel('x'),ylabel('y'),zlabel('z')
axis equal
%% Form the sheets for the location
filename = "Mic_positions.txt"
lines = strings(2*length(thetas)+1,1);
for i = 1:length(thetas)
    % alinged with karaub
    minor_loc = z_plane_locs(i,:);
    major_loc = y_plane_locs(i,:);
    
    Minor=sprintf(name_fmt,'Minor',i,minor_loc(1),minor_loc(2)...
        ,minor_loc(3),minor_loc(4));
    Major=sprintf(name_fmt,'Major',i,major_loc(1),major_loc(2)...
        ,major_loc(3),major_loc(4));
    lines(i) = Minor;
    lines(length(thetas)+i+1) = Major;
end
writelines(lines,filename)
%% Block size guide and strouhal number calculations
fs = 204.8E3;%Hz
T = 5; 
N_tot = fs*T
df = 50; %Hz
N = fs/df
f_max = fs/2 - df %Hz
St_fmax = f_max * De/Uj

fs_s1 = fs * h/c_inf
%we are sampling way quicker than what is needed to match
fs_sc = fs * De/Uj * Ma/De_s

dts = 1/fs_sc
step = dts/1e-3
%let the sampling frquency match

%%

dts = 0.07; % 0.05 %this is the fwh-sampling interval, not simulation interval
fs_s = 1/dts
St_s = fs_s * De_s/Ma %Ma since c_inf_sim = 1 

fs_p = St_s * Uj/De

fs_p1= fs_s * c_inf/h
St_smax = St_s/2 - St_bw

dfs_bw = St_bw*Ma/De_s
blocksize_s = fs_s/dfs_bw


t_block_s = blocksize_s*dts
%need 67.3650 to generate one block of ~4000 samples i.e. 3 days
%very likely to not have a lot of samples to do things... (I am a bit
%concerned...)

%%
St_scr = 0.37; %estimate of screech frequency
fs_scr = St_scr * Ma/De_s
l_scr = 1/fs_scr
delta_scr = De_s/(8*St_scr*Ma)

St_cutoff_high = 4;
fs_cutoff_high = St_cutoff_high*Ma/De_s;
delta_high = De_s/(8*St_cutoff_high*Ma)

%% Determining the minimum level necessary for capturing screech
fprintf("upstream spacing determination\n")
HCP_delta = 2.56
LEVEL = log2(HCP_delta/((c_inf*De)/(8*St_scr*Uj*h)))
%% 
%St_tgt = 0.1 ; %rough range where information is currently missing from
%fs_tgt = St_tgt * Ma/De_s
%ls_tgt = 1/fs_tgt