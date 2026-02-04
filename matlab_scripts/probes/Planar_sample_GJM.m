%% k_gjm based estimation for sampling
clear;clc;close all;
%% compute the jet Mach number and jet velocity
NPR = 3;
gamma = 1.4;
R = 1/gamma;
Ta = 1;
c_inf = sqrt(gamma*R*Ta)
%we know that c_inf = 1
%therefore given c = sqrt(gamma R T) and normailzed rho, R and T are all one
%
Mj = sqrt(2/(gamma-1)*(NPR^((gamma-1)/gamma)-1)) %jet Mach number
Tj = (1+ (gamma-1)/2*Mj^2)^(-1)
cj = sqrt(gamma*R*Tj)
Ma = Mj*cj/c_inf %acoustic Mach number

%% from measurements we have
lambda_shock = 1.3; %nondim wavelength of shockcell structure from paraview measurement
k_shock = 2*pi/(lambda_shock)


c_kh = 0.8*Ma %nondim phase velocity of Kelvin Helmholtz wave

L = 2;
h = 1;
De = 1.3*(L*h)^(0.625)/(L+h)^(0.25)

St_fact = De/Ma; %factor to go from nondim freq to Strouhaul
St_screech = 0.37; %from post-processing
f_screech = St_screech/St_fact
omega_screech = 2*pi*f_screech
k_screech = omega_screech/c_inf
lambda_screech = 2*pi/k_screech
dx_screech = lambda_screech/16

k_kh = (2*pi*f_screech)/(c_kh)

k_gjm = k_kh - k_shock
omega_gjm = c_inf/k_gjm
f_gjm = omega_gjm/2*pi
lambda_gjm = 2*pi/k_gjm
dx_gjm = lambda_gjm/16

%% Sampling time adjustment
T_screech = 1/f_screech
dt_screech = T_screech/22
dt_fwh_samp = 50
% finalized dx and dt are
deltax = floor(abs(dx_gjm)/0.1)*0.1
deltat = floor(abs(dt_screech)/0.01)*0.01
dt_samp = round(deltat*1000/dt_fwh_samp)*dt_fwh_samp
%% build the nozzle interior grid 

%% xz plane interior
%this is the CD nozzle section
x_min = 3.1; %slightly upstream on the end of the converging section
x_tr1 = -2.61052632; %start of converging section
x_tr2 = -1.0152632; %start of diverging section
x_tr3 = 0; %outlet
x_noz = -x_min +mod(x_min,deltax):deltax:0;
Nx_noz = length(x_noz);
%build the necessary profile
ind1 = find(x_noz - x_tr1 > 0,1)
ind2 = find(x_noz - x_tr2 > 0,1)
%ind3 -> x = 0, is the end of x_noz
%% determining the height of the section
z_1 = 1; %half height of entry to converging section
z_2 = 0.42498737; %half height of throat
z_3 = 0.5; %half height of outlet

%introducing the one cell shift
dl_base = 3.2
dl_z    = dl_base/(2^8)

%shrinking the size
z_1 = z_1 - dl_z
z_2 = z_2 - dl_z
z_3 = z_3 - dl_z
%% 

% I need to build the bounds
z_profile_top = ones(1,Nx_noz);
slope1 = (z_2 - z_1)/(x_tr2 - x_tr1);
z_profile_top(ind1:ind2-1) = (x_noz(ind1:ind2-1)-x_tr1)*slope1 + z_1;
slope2 = (z_3 - z_2)/(x_tr3 - x_tr2);
z_profile_top(ind2:end) = (x_noz(ind2:end) - x_tr2)*slope2 + z_2;
z_profile_bot = -z_profile_top;
%validate profile
figure; hold on; plot(x_noz, z_profile_top,'ro--'); axis equal
plot(x_noz,z_profile_bot,'bs--')   

%I have top and bottom boundaries, now grid the system such that the grid
%always satisfyies the condition;
z_noz_init = -0.5:deltax:0.5; %using outlet as bound
Nz_noz = length(z_noz_init);
xz_nozzle = zeros((Nx_noz)*Nz_noz,3);
for i = 1:Nx_noz
    z_min = z_profile_bot(i);
    z_max = z_profile_top(i);
    z_cur = linspace(z_min,z_max,Nz_noz);
    for k = 1:Nz_noz
        ind = (i-1)*Nz_noz + k;
        xz_nozzle(ind,1) = x_noz(i);
        xz_nozzle(ind,3) = z_cur(k);
    end
end

X_noz = reshape(xz_nozzle(:,1),[Nx_noz, Nz_noz]);
Z_noz = reshape(xz_nozzle(:,3),[Nx_noz, Nz_noz]);
%validation check
figure; hold on;
scatter(xz_nozzle(:,1),xz_nozzle(:,3))

%% xy plane interior

y_noz = (-1+dl_z):deltax:0;
y_noz = [y_noz,0,flip(abs(y_noz))];
Ny_noz = length(y_noz);
xy_nozzle = zeros((Nx_noz)*Ny_noz,3);

for i = 1:Nx_noz
    for j = 1:Ny_noz
        ind = (i-1)*Ny_noz + j;
        xy_nozzle(ind,1) = x_noz(i);
        xy_nozzle(ind,2) = y_noz(j);
    end
end
%% external planes external
%for x = 0:85; // Do I need to collect this long? probably not!
%for y = +/-7 or z = +/-7
x_start = deltax;
%x_end = 85;
x_end = 24; %just beyond the end of the potential core
y_start = -7;
y_end = 7;
x = x_start:deltax:x_end;
y = y_start:deltax:y_end;
z_start= y_start;
z_end = y_end;
z = z_start:deltax:z_end;

Nx = length(x);
Ny = length(y);
Nz = length(z);

xy_grid = zeros(Nx*Ny,3);
xz_grid = zeros(Nx*Nz,3);
for i = 1:Nx
    for j = 1:Ny
        ind = (i-1)*Ny + j;
        xy_grid(ind,1) = x(i);
        xy_grid(ind,2) = y(j);
    end
end

for i = 1:Nx
    for k = 1:Nz
        ind = (i-1)*Nz + k;
        xz_grid(ind,1) = x(i);
        xz_grid(ind,3) = z(k);
    end
end

%% pad additional points external to the geometry

%% xz plane
%determine number of points needed outside
Nz_ex = Nz - Nz_noz %excess number of points
Nz_ex_half= Nz_ex/2
%building the z_profile as before (x_min = 3.1)
x_tr1 = 0;
x_tr2 = -1.412;
x_tr3 = -1.899;

%top profile with additional add on
z_tr1 = 0.584 + dl_z;
z_tr2 = 1.501 + dl_z;
z_tr3 = 1.676 + dl_z;

%build z_bot profile
z_top = z_end;
x_noz_f = flip(x_noz);
ind2 = find(x_noz_f - x_tr2 < 0 ,1)
ind3 = find(x_noz_f - x_tr3 < 0 ,1)
slope1 = (z_tr2 - z_tr1)/(x_tr2 - x_tr1)
slope2 = (z_tr3 - z_tr2)/(x_tr3 - x_tr2)
z_bot = zeros(1,Nx_noz);
z_bot(1:ind2-1) = z_tr1 + slope1*(x_noz_f(1:ind2-1) - x_tr1)
z_bot(ind2:ind3-1) = z_tr2 + slope2*(x_noz_f(ind2:ind3-1) - x_tr2)
z_bot(ind3:end) = z_tr3

z_bot = flip(z_bot);
figure
plot(x_noz, z_bot, 'r-o')

%build the grid as needed
xz_extern = zeros(Nx_noz*Nz_ex,3); 
for i = 1:Nx_noz
    z_profile = linspace(z_bot(i),z_top,Nz_ex_half);
    z_profile = [-flip(z_profile),z_profile]; %get the negative side as well
    for k = 1:Nz_ex
        ind = (i-1)*Nz_ex + k;
        xz_extern(ind,1) = x_noz(i);
        xz_extern(ind,3) = z_profile(k);
    end
end


%% xy plane
Ny_ex = Ny - Ny_noz
Ny_ex_half = Ny_ex/2

x_tr1 = -3.179;
x_tr2 = 0; 

y_tr1 = 1.347+dl_z;
y_tr2 = 1.074 + dl_z;

slope = (y_tr2 - y_tr1)/(x_tr2 - x_tr1)
y_bot = slope*(x_noz-x_tr1)+y_tr1;
y_top = y_end;
figure;
plot(x_noz,y_bot,'r-o')
% build the grid as needed
xy_extern = zeros(Nx_noz*Ny_ex,3);
for i = 1:Nx_noz
    y_profile = linspace(y_bot(i), y_top,Ny_ex_half);
    y_profile = [-flip(y_profile),y_profile];
    for j = 1:Ny_ex
        ind = (i-1)*Ny_ex + j;
        xy_extern(ind,1) = x_noz(i);
        xy_extern(ind,2) = y_profile(j);
    end
end
%%
fv = stlread("nozzle_geometry.stl");
vertices = fv.Points;
factor = 1/0.475;
scaledVertices = vertices*factor;


%% validation with nozzle geometry
figure; hold on; view(2)
scatter3(xy_grid(:,1),xy_grid(:,2),xy_grid(:,3),'r','filled')

%scatter3(xz_grid(:,1),xz_grid(:,2),xz_grid(:,3),'b','filled')
scatter3(xy_nozzle(:,1),xy_nozzle(:,2),xy_nozzle(:,3),'m','filled')
scatter3(xy_extern(:,1),xy_extern(:,2),xy_extern(:,3),'m','filled')
%scatter3(xz_nozzle(:,1),xz_nozzle(:,2),xz_nozzle(:,3),'g','filled')
%scatter3(xz_extern(:,1),xz_extern(:,2),xz_extern(:,3),'g','filled')
trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), scaledVertices(:, 3),'FaceColor','c','FaceAlpha',0.1,'EdgeColor','none')
axis equal; 
% 
figure; hold on; view(0,0)
%scatter3(xy_grid(:,1),xy_grid(:,2),xy_grid(:,3),'r','filled')

scatter3(xz_grid(:,1),xz_grid(:,2),xz_grid(:,3),'b','filled')
%scatter3(xy_nozzle(:,1),xy_nozzle(:,2),xy_nozzle(:,3),'m','filled')
%scatter3(xy_extern(:,1),xy_extern(:,2),xy_extern(:,3),'m','filled')
scatter3(xz_nozzle(:,1),xz_nozzle(:,2),xz_nozzle(:,3),'g','filled')
scatter3(xz_extern(:,1),xz_extern(:,2),xz_extern(:,3),'g','filled')
trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), scaledVertices(:, 3),'FaceColor','c','FaceAlpha',0.1,'EdgeColor','none')
ylim([-8,8])
axis equal; 

%% append of data

xy_upstream_combined = zeros(Nx_noz*Ny,3);

for i = 1:Nx_noz
    %for each x position
    ind_str = 1 + Ny*(i-1);
    ind_sp1 = ind_str + Ny_ex_half;
    ind_sp2 = ind_sp1 + Ny_noz;
    ind_end = Ny*i;
    nozzle_row_inds = ((i-1)*Ny_noz+1):i*Ny_noz;
    extern_low_inds = (((i-1)*Ny_ex)+1):((i-1)*Ny_ex+Ny_ex_half);
    extern_hi_inds  = ((i-1)*Ny_ex+Ny_ex_half+1):(i*Ny_ex);
    xy_upstream_combined(ind_str:(ind_sp1-1),:) = xy_extern(extern_low_inds, :);
    xy_upstream_combined(ind_sp1:(ind_sp2-1),:) = xy_nozzle(nozzle_row_inds,:);
    xy_upstream_combined(ind_sp2:ind_end,:) = xy_extern(extern_hi_inds,:);
end

xy_grid = [xy_upstream_combined;xy_grid];


xz_upstream_combined = zeros(Nx_noz*Nz,3);

for i = 1:Nx_noz
    %for each x position
    ind_str = 1 + Nz*(i-1);
    ind_sp1 = ind_str + Nz_ex_half;
    ind_sp2 = ind_sp1 + Nz_noz;
    ind_end = Nz*i;
    nozzle_row_inds = ((i-1)*Nz_noz+1):i*Nz_noz;
    extern_low_inds = (((i-1)*Nz_ex)+1):((i-1)*Nz_ex+Nz_ex_half);
    extern_hi_inds  = ((i-1)*Nz_ex+Nz_ex_half+1):(i*Nz_ex);
    xz_upstream_combined(ind_str:(ind_sp1-1),:) = xz_extern(extern_low_inds, :);
    xz_upstream_combined(ind_sp1:(ind_sp2-1),:) = xz_nozzle(nozzle_row_inds,:);
    xz_upstream_combined(ind_sp2:ind_end,:) = xz_extern(extern_hi_inds,:);
end

xz_grid = [xz_upstream_combined;xz_grid];
figure; hold on; view(0,0)
scatter3(xz_grid(:,1),xz_grid(:,2),xz_grid(:,3),'b','filled')
trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), scaledVertices(:, 3),'FaceColor','c','FaceAlpha',0.1,'EdgeColor','none')
ylim([-8,8])
axis equal; 



%% injection grid
inject_grid = xz_grid;
inject_grid(:,2) = 0.143; %first injection port center;
figure; hold on; view(0,0)
scatter3(inject_grid(:,1),inject_grid(:,2),inject_grid(:,3),'b','filled')
trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), scaledVertices(:, 3),'FaceColor','c','FaceAlpha',0.1,'EdgeColor','none')
ylim([-8,8])
axis equal; 




%% export 
% writematrix(xy_grid,'xy_plane_gjm_noheader.txt')
% writematrix(xz_grid,'xz_plane_gjm_noheader.txt')
% writematrix(inject_grid,'inject_plane_gjm_noheader.txt')
% writematrix(xy_nozzle,'xy_plane_nozzle_gjm_noheader.txt')
% writematrix(xz_nozzle,"xz_plane_nozzle_gjm_noheader.txt")
% 



