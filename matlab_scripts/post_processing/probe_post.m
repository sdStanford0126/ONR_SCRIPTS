%axial probe initail processing
%following Olivia Martin's implementation
%done by Steven Dai
clear;clc;close all
%% setup input parameters
input_dir = "./pcprobes/";
output_dir = "./";
start_id = 309550;
end_id = 423000;
dt = 50;
tid = start_id:dt:end_id;

prefix = "axial.";
suffix = ".pcd";
dataname_fmt = prefix + "%08d" + suffix;

pts_file = input_dir +prefix+"pxyz"
%% grab position data
xyz = importdata(pts_file).data;
x = xyz(:,2);
y = xyz(:,3);
z = xyz(:,4);
% only extract the centerline
ind_center = find(xyz(:,3)==0 & xyz(:,4)==0);
xc = xyz(ind_center,2); %checked that it is indeed inline
yc = xyz(ind_center,3);
zc = xyz(ind_center,4);
%% grab centerline pressure data
%pressure correspond to fourth column
% p_hist = zeros(length(tid),length(ind_center));
% for i = 1:length(tid)
%     tid_cur = tid(i);
%     dataname = input_dir+sprintf(dataname_fmt,tid_cur);
%     probedata = importdata(dataname).data;
%     p_hist(i,:) = probedata(ind_center,4);
%     %i
% end
%% save/load data for future use
%save("axial_centerline_p.mat","p_hist","xc")
load("axial_centerline_p.mat")
%% need to get
f2St = 1.3126; %upscaling factor
t_step = 0.001;
dt = 50;
deltaT = t_step*dt;
fs = 1/deltaT;
N = length(tid);
freqRange = fs/N*(0:N-1);
StRange = freqRange*f2St; %rescale to St
p_fft = 1/N*fft(p_hist,[],1); %take fft in temporal direction
P_PSD = 2.*abs(p_fft).^2;

%%
figure;
%loglog(freqRange(1:end/2),2*abs(p_spect(1:end/2,1)),DisplayName=sprintf("x=%2d",xc(1)))

for i = 2:10:length(xc)
   loglog(freqRange(1:end/2),P_PSD(1:end/2,i),DisplayName=sprintf("x=%2d",xc(i)),LineWidth=1)
   hold on;
end
legend
xlim([0.05,0.3])
