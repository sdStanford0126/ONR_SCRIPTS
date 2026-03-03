%visualize init SPOD results
%Steven Dai, based on examples provided by O. T. Schmidt (oschmidt@ucsd.edu)
clear;clc;close all; setPlotPref(1.5,'latex',14)
%% laod data SPOD data
addpath("baseline_new/130M/")

%% load xy plane (minor observer plane)
input_dir = "~/Downloads/spod_results_xz_p/nfft390_novlp193_nblks6/";
data_fname = "spod_energy.mat";
plane_name = "Major_observer";
SPOD_eng_plot_name = "SPOD_p_var_" + plane_name + ".png"

load(input_dir + data_fname);

%% some basline configurations
rho = 1.2; %kg/m^3
c   = 346.1; %m/s %sufficiently close (346.0613
Punit = rho*c^2;

NPR = 3;
gamma = 1.4;
Ma = jetMach_isen_ideal(NPR,gamma)%jet mach number acoustic 
h_s = 1;
L_s = 2*h_s;
De_s = 1.3*(L_s*h_s)^(0.625)/(L_s+h_s)^(0.25)
fact_St = De_s/Ma

De = De_s * 0.475 /39.37; %m
h = De/De_s;
fact_St_p = De / (Ma*c)
%% plotting
St = f*fact_St;
figure;
loglog(St,L)
xlabel("$St$")
ylabel("$|p'|^2/St$") 
title("SPOD P VAR, " + plane_name)
saveas(gcf,SPOD_eng_plot_name)

%% direct visualization of the data
mode_fmt = "spod_f%04d.mat";

%need to find the screech frequency number
df = f(2)-f(1);
St_sc = 0.358;
%St_sc = 0.1;
%St_sc = 0.1
f_sc = St_sc/fact_St
fi = find(abs(f-f_sc)<df/2) %index of frequency

%load appropriate mode shape
mode_fname = sprintf(mode_fmt, fi)
load(input_dir+mode_fname)
%loaded variable is Psi
del_x = 0.1;
X = repmat(transpose(0:del_x:85),[1,141]);
Z = repmat(-7:del_x:7,[851,1]);
%do the same visualization of this mode
nt = 30;
Ps1 = Psi(:,:,1); %take first mode
draw_first_mode(fi,f,Ps1,nt,X,Z,fact_St,L)


%% lets do streamwise decomposition
%given what we know we should get the mode of interest and decompose
Len = 85; %extent in x/h
Nx = 851;
kx = transpose(2*pi/Len * (-(Nx-1)/2:(Nx-1)/2)); %the wavenumbers
z_vec = -7:del_x:7;
Kx = repmat(kx,[1,141]);
mi = 1;
Ps1 = Psi(:,:,mi); %take only the first mode
%windowing of data
w_hann = repmat(hann(Nx),[1,141]);
figure
plot(hann(Nx))
Ps1_hat = fft(w_hann.*Ps1,[],1); %take fft in the x direction, note the windowing

Ps1_hat = 2*fftshift(Ps1_hat,1); %do fft shift as needed

k_acoustics = [-(2*pi*f_sc),(2*pi*f_sc)]


figure;
contour(-Kx,Z,abs(Ps1_hat),50)
xline(k_acoustics,'r-',LineWidth=1.5)
xlim([-5,5]), ylim([-4.5,4.5]),xlabel("$k_x h$"), ylabel("$y/h$")
colorbar()
saveas(gcf,'screech_mode_modulus.png')
%% now isolate the upstream/downstream traveling waves
ind_sep = find(kx==0)
%% upstream traveling (due to inversion of frequencies)

Ps1_hat_p = Ps1_hat;
Ps1_hat_p(1:ind_sep,:) = 0; %zero out upstream traveling waves;
Ps1_p = ifft(ifftshift(Ps1_hat_p,1),[],1); % reform into physical space
draw_first_mode(fi,f,Ps1_p,30,X,Z,fact_St,L)
figure;
pcolor(X,Z,real(Ps1_p)), shading interp, axis equal tight
xlim([0,15]), ylim([-4.5,4.5]),xlabel("$x/h$"), ylabel("$y/h$")
saveas(gcf,"upstream_mode_shape.png")

%% downstream traveling (due to inversion of frequencies)
Ps1_hat_m = Ps1_hat;
Ps1_hat_m(ind_sep:end,:) = 0; %zero out downstream traveling waves;
Ps1_m = ifft(ifftshift(Ps1_hat_m,1),[],1); % reform into physical space
draw_first_mode(fi,f,Ps1_m,nt,X,Z,fact_St,L)

figure;
pcolor(X,Z,real(Ps1_m)), shading interp, axis equal tight
xlim([0,15]), ylim([-4.5,4.5]),xlabel("$x/h$"), ylabel("$y/h$")
saveas(gcf,"downstream_mode_shape.png")
%% helper function
function draw_first_mode(fi,f,Ps1,nt,X,Z,fact_St,L)
T = 1/f(fi);
time = linspace(0,3*T,nt);
figure;

for ti = 1:nt
    pcolor(X,Z,real(Ps1*exp(2i*pi*f(fi)*time(ti)))), shading interp, axis equal tight, caxis(max(abs(caxis))*[-1 1])
            xlabel('x'), ylabel('r'), 
            title(['St=' num2str(f(fi)*fact_St,'%.2f') ', mode ' num2str(1) ', ' ...
                '$\lambda$=' num2str(L(fi,1),'%.2g')]) %only the first mode
            xlim([0,25]), ylim([-7,7])
    drawnow
    hold off
end
end
