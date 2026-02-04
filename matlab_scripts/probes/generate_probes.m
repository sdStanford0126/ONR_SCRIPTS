% Write probe points to data file (from Olivia, modified by Steven)
clear all
format long %add to all other probe generation file for future use

fname = strcat('axial_profiles_probes.txt');
addpath("../")
%setPlotPref(1.4,'latex',14)

%% quick pre-amble
x_pl = [0,0.425,1.84768,3.55153,6.47314,12.6154, 19.0734, 25, 31, 40];
y_pl = [1,1.01648,1.1,1.17,1.29,1.85,2.15,2.85,3,3.8];
z_pl = [0.5,0.50541,0.657779, 0.838731,1.27942,2.2553,2.92125,3.50369,3.58816,3.75443];
x_samp = [0,0.5,1.8,3.5,6.5,13,19,24,31,40]; %x sampling positions 
figure; hold on; box on;
plot(x_pl,y_pl,'r-o',DisplayName="Major axis plane")
plot(x_pl,z_pl,'b-^',DisplayName="Minor axis plane")
legend(location="northwest",AutoUpdate="off")
xline(x_samp,'k--',DisplayName=" ")
xlabel("x/h")
ylabel("y/h or z/h")

saveas(gcf,'Jet_axis_switch_baseline_144M.png')
%% mesh spacings
%dx = 0.1;
dy = 0.05;
dz = 0.05;

%x = [0:dx:40];
x = [0,0.5,1.8,3.5,6.5,13,19,24,31,40]; %x sampling positions 
y = [-10:dy:10];
z = [-10:dz:10];

Nx = length(x);
Ny = length(y);
Nz = length(z);

%% create probes

dataMat = zeros(Nx*Ny*Nz, 3);
[X,Y,Z] = meshgrid(x,y,z);

dataMat(:,1) = X(:);
dataMat(:,2) = Y(:);
dataMat(:,3) = Z(:);
% 
% fid = fopen(fname,'w');
% fprintf(fid,'X\tY\tZ\n');
% for k=1:size(dataMat, 1)
%    fprintf(fid,'%.8f\t%.8f\t%.8f\n',dataMat(k,:));
% end

%fclose(fid);


%%
figure
plot(dataMat(:,1))
figure
plot(dataMat(:,2))
figure
plot(dataMat(:,3))
