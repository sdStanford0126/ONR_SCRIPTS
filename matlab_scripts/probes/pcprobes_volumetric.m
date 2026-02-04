%Volumetric Sampling
%Could I do this with one timestep? 
%worth a shot
clear;clc;close all;
outdir = "/Users/steven/Downloads/";
%% sampling mesh choice
HCP_DELTA=2.56
LEVEL = 9 %target level 6 
deltax = HCP_DELTA/(2^LEVEL)
%% nozzle interior
%this is the CD nozzle section, up to the throat
x_noz = -1.0105263 +mod(1.0105263,deltax):deltax:0;
Nx_noz = length(x_noz);
%build the necessary profile
%ind1 = find(x_noz + 2.61052632 > 0,1)
ind2 = find(x_noz + 1.0152632 > 0,1)
%x_tr1 = -2.61052632; %start of converging section
x_tr2 = -1.0152632; %start of diverging section
x_tr3 = 0; %outlet
%z_1 = 1; %half height of entry to converging section
z_2 = 0.42498737-deltax/2; %half height of throat
z_3 = 0.5-deltax/2; %half height of outlet
% I need to build the bounds
slope2 = (z_3 - z_2)/(x_tr3 - x_tr2);
z_profile_top = (x_noz - x_tr2)*slope2 + z_2;
z_profile_bot = -z_profile_top;
%validate profile
figure; hold on; plot(x_noz, z_profile_top,'ro--'); axis equal
z_noz_init = -0.5+deltax:deltax:0.5-deltax; %this is upstream of the throat
Nz_noz = length(z_noz_init);
y_noz = -1+deltax:deltax:1-deltax;
Ny_noz = length(y_noz);
xyz_nozzle = zeros(Nx_noz*Ny_noz*Nz_noz,3);
for i = 1:Nx_noz
    z_min = z_profile_bot(i);
    z_max = z_profile_top(i);
    z_cur = linspace(z_min,z_max,Nz_noz);
    for j = 1:Ny_noz
        for k = 1:Nz_noz
            ind = (i-1)*Nz_noz*Ny_noz + (j-1)*Nz_noz + k;
            xyz_nozzle(ind,1) = x_noz(i);
            xyz_nozzle(ind,3) = z_cur(k);
            xyz_nozzle(ind,2) = y_noz(j);
        end
    end
end

X_noz = reshape(xyz_nozzle(:,1),[Nx_noz,Ny_noz ,Nz_noz]);
Y_noz = reshape(xyz_nozzle(:,2), [Nx_noz,Ny_noz ,Nz_noz]);
Z_noz = reshape(xyz_nozzle(:,3),[Nx_noz, Ny_noz ,Nz_noz]);
%% validation check
figure;
fv = stlread("nozzle_geometry.stl");
vertices = fv.Points;
factor = 1/0.475;
scaledVertices = vertices*factor; 
hold on;
trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), scaledVertices(:, 3),'FaceColor','r','FaceAlpha',0.3)
scatter3(xyz_nozzle(:,1),xyz_nozzle(:,2),xyz_nozzle(:,3))
axis equal
%%
%writematrix(xyz_nozzle,outdir+"nozzle_interior_vol_L9_no_header.txt",Delimiter=" ")
%% shockcell region downstream
% from paraview, established that we would need at least x = 15 to capture,
% extend to x = 20 
% with extend +/- 2
% volume would be 2 * 2 *20 = 80 h^3 ==> giving 20 million points...

% all shockcells 

y = -2:deltax:2;
z = -2:deltax:2;
x = deltax:deltax:30+deltax;
mesh_size = length(y)*length(x)*length(z);
Ny = length(y);
Nx = length(x);
Nz = length(z);



%% Multiple plane injection port sampling
clear("xyz_nozzle")

BL = 2.56;
LEVEL = 11;

Dx = BL/2^(LEVEL)
x_shift = 0.5*Dx;
x_str = -0.28027368+x_shift
x_end = -0.22498947-x_shift
dy = 0.09212632 - 2*x_shift; %width for each port
ys_str_pos = [0.0968+x_shift,0.38250526+x_shift,0.66820421+x_shift];
ys_str = [flip(-(ys_str_pos+dy)),ys_str_pos]
numPort = length(ys_str);

z_shift = 2*Dx; %two cell inward
z_bot_str = -0.47919493 - z_shift;
z_bot_end = -0.48329875 - z_shift;

Nx = round((sqrt((x_end-x_str)^2 + (z_bot_str - z_bot_end)^2))/Dx)
Nz = Nx
Ny = round(dy/Dx)
xs = linspace(x_str,x_end,Nx);



zs_bot= linspace(z_bot_str,z_bot_end,Nz);
zs_top = -zs_bot;

Nlayers = 4;

dz = 0.015;
zs_bot = repmat(zs_bot,Nlayers,1);
zs_top = repmat(zs_top,Nlayers,1);

for i = 2:Nlayers
    zs_bot(i,:) = zs_bot(i,:) - dz*(i-1);
    zs_top(i,:) = zs_top(i,:) + dz*(i-1);
end

%% Top port
%
pos_top_ports = zeros(Nx*Ny*Nlayers*numPort,3); %all top port locations
%arrangement from -y to +y
for i = 1:numPort
    y_str = ys_str(i);
    ys = linspace(y_str,y_str+dy,Ny);
    for j = 1:Nlayers
        for m = 1:Ny
            for n = 1:Nx
                %for each individual port, go from layer closest to the
                %diverging plane to plenum
                ind = (i-1) * (Nx*Ny*Nlayers) + (j-1)*(Nx*Ny) ...
                    + (m-1)*Nx + n; %%local index
                x = xs(n);
                y = ys(m);
                z = zs_top(j,n); % Nlayer,Nx
                pos_top_ports(ind,:)=[x,y,z];
            end
        end
    end
end
%validate
figure;hold on;
trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), ...
    scaledVertices(:, 3),'FaceColor','r','FaceAlpha',0.3)
scatter3(pos_top_ports(:,1),pos_top_ports(:,2),pos_top_ports(:,3))
axis equal;

writematrix(pos_top_ports,outdir+"ports_layered_top_no_header.txt")
%% Bot port
%
pos_bot_ports = zeros(Nx*Ny*Nlayers*numPort,3); %all top port locations
%arrangement from -y to +y
for i = 1:numPort
    y_str = ys_str(i);
    ys = linspace(y_str,y_str+dy,Ny);
    for j = 1:Nlayers
        for m = 1:Ny
            for n = 1:Nx
                %for each individual port, go from layer closest to the
                %diverging plane to plenum
                ind = (i-1) * (Nx*Ny*Nlayers) + (j-1)*(Nx*Ny) ...
                    + (m-1)*Nx + n; %%local index
                x = xs(n);
                y = ys(m);
                z = zs_bot(j,n); % Nlayer,Nx
                pos_bot_ports(ind,:)=[x,y,z];
            end
        end
    end
end
%validate
figure;hold on;
trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), ...
    scaledVertices(:, 3),'FaceColor','r','FaceAlpha',0.3)
scatter3(pos_bot_ports(:,1),pos_bot_ports(:,2),pos_bot_ports(:,3))
axis equal;
view(3)

writematrix(pos_bot_ports,outdir+"ports_layered_bot_no_header.txt")

%% volumetric sampling around two central ports at bottom
clear("pos_top_ports","pos_bot_ports")

BL = 2.56
LEVEL = 10 %this is already an over estimate

Dx = BL/2^(LEVEL); %targe resolution
deltax = Dx*0.5;
%x_bounds
x_str = -0.35;
x_end = 1;
x_in = x_str:Dx:0;
Nx_in = length(x_in);
x_out = Dx:Dx:x_end;
Nx_out = length(x_out);
xs = [x_in,x_out];
Nx = length(xs);

%y_bounds
y_str = 0;
y_end = 0.475+0.194/2;
ys = y_str:Dx:y_end;
Ny = length(ys);
%need to deal with the z transition
z_top = 0;
z_bot_out = -0.5+deltax;
z_bot_str = -(slope2*(x_str-x_tr2)+0.425-deltax);
z_bot_in = linspace(z_bot_str,z_bot_out,Nx_in);
z_bot = [z_bot_in,ones(1,Nx_out)*z_bot_out];
zs_out = z_bot_out:Dx:z_top;
Nz = length(zs_out);

probe_inject_port = zeros(Nx*Ny*Nz,3);
%move from 0 to +y
for i = 1:Nx
    %for a given x location generate correct z
    zs = linspace(z_bot(i),z_top,Nz);
    for j = 1:Ny
        for k = 1:Nz
            ind = (j-1)*Nx*Nz + (k-1)*Nx + i;
            probe_inject_port(ind,:) = [xs(i),ys(j),zs(k)];
        end
    end
end

%% 
figure;hold on;
trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), ...
    scaledVertices(:, 3),'FaceColor','r','FaceAlpha',0.3)
scatter3(probe_inject_port(:,1),probe_inject_port(:,2),probe_inject_port(:,3))
axis equal;
view(2)

%%
writematrix(probe_inject_port,outdir+"injection_port_nearfield_no_header.txt")
%% multiple plane sampling
y_sample_loc = -1:0.1:1;
z_smaple_loc = -0.5:0.1:0.5;

Ny_samp = length(y_sample_loc);
Nz_samp = length(z_smaple_loc);
%setup xy_plane sampleing
xy_sc_plane = zeros(Nx*Ny*Nz_samp,3); 
for k = 1:Nz_samp
    for j = 1:Ny
        for i = 1:Nx
            ind = (k-1)*Nx*Ny + (i-1)*Ny + j;
            xy_sc_plane(ind,1) = x(i);
            xy_sc_plane(ind,2) = y(j);
            xy_sc_plane(ind,3) = z_smaple_loc(k);
        end
    end
end

figure; scatter3(xy_sc_plane(:,1),xy_sc_plane(:,2),xy_sc_plane(:,3));
xlabel('x');ylabel('y');zlabel('z'); axis equal;

%setup xz_plane sampling
xz_sc_plane = zeros(Nx*Ny_samp*Nz,3);
for j = 1:Ny_samp
    for k = 1:Nz
        for i = 1:Nx
            ind = (j-1)*Nx*Ny + (i-1)*Nz + k;
            xz_sc_plane(ind,1) = x(i);
            xz_sc_plane(ind,2) = y_sample_loc(j);
            xz_sc_plane(ind,3) = z(k);
        end
    end
end
figure; scatter3(xz_sc_plane(:,1),xz_sc_plane(:,2),xz_sc_plane(:,3));
xlabel('x'); ylabel('y'); zlabel('z'); axis equal;

writematrix(xy_sc_plane,'shock_xy_planes_no_header.txt', Delimiter=" ")
writematrix(xz_sc_plane,'shock_xz_planes_no_header.txt', Delimiter=" ")
% xyz = zeros(Nx*Ny*Nz,3);
% 
% for i = 1:Nx
%     for j = 1:Ny
%         for k = 1:Nz
%             ind = (i-1)*Ny*Nz + (j-1)*Nz + k;
%             xyz(ind,1) = x(i);
%             xyz(ind,2) = y(j);
%             xyz(ind,3) = z(k); 
%         end
%     end
% end
% 
% 
% figure;
% scatter3(xyz(:,1),xyz(:,2),xyz(:,3))
% axis equal
% 
% % 
% 
% writematrix(xyz,'shockcell_vol_fine_no_header.txt',Delimiter=" ")




