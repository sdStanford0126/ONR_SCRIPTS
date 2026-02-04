%% probes_location
% Date: Nov. 30th, 2024
clear;clc;close all
%% setup probe locations
%separate files for separate target zones to handle post-processing

%% TODO: probe locations within the main nozzle divergent section
% from x/h = -1.01052632 to x/h = 0

%define start and end position for x and z
x_str = -1.01052632;
x_end = 0.0;

%need to introduce a z_shift to avoid
%typically LEVEL 7 should suffice
z_shift = 2.56/(2^9);

z_str_bot = -0.42498737 + z_shift;
z_str_top = 0.42498737 - z_shift;

z_end_bot = -0.5 + z_shift;
z_end_top = 0.5 - z_shift; 

y_lim = 1.0;

%find the slope 
z_slope_top = (z_end_top - z_str_top)/(x_end - x_str);
z_slope_bot = (z_end_bot - z_str_bot)/(x_end - x_str);

%setting equal position to make later analysis easier (i.e. uniform
delta = 2.56/(2^7);

xs = x_str:delta:x_end;
Nx = length(xs);
%y position declared ahead since no change in behavior
ys = -y_lim:delta:y_lim;
Ny = length(ys);
%z position extends into nozzle (for now)
zs = z_end_bot:delta:z_end_top; 
Nz = length(zs);
Positions = zeros(Nx*Ny*Nz,3);

%%
figure
scatter3(Positions(:,1),Positions(:,2),Positions(:,3))
%%
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            x = xs(i);
            y = ys(j);
            z = zs(k);
            ind = (i-1)*Ny*Nz + (j-1)*Nz + k;
            Positions(ind,1) = x;
            Positions(ind,2) = y;
            Positions(ind,3) = z;
        end
    end
end
%% read nozzle geometry for plotting
fv = stlread("nozzle_geometry.stl");
vertices = fv.Points;
factor = 1/0.475;
scaledVertices = vertices*factor;

%% 
figure;
scatter3(Positions(:,1),Positions(:,2),Positions(:,3)); hold on;
trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), scaledVertices(:, 3),'FaceColor','r','FaceAlpha',0.3)
daspect([1 1 1])
%% 
writematrix(Positions,'interior_probes_xyz_no_header.txt',Delimiter=" ")

%% TODO: probe locations at the injection ports (top and bottom)
x_shift = 0.5*3.2/2^11;
x_str = -0.28027368+x_shift
x_end = -0.22498947-x_shift
dy = 0.09212632 - 2*x_shift;
ys_str_pos = [0.0968+x_shift,0.38250526+x_shift,0.66820421+x_shift];
ys_str = [flip(-(ys_str_pos+dy)),ys_str_pos]


% 200 points per port, 12 ports give 2400 pts 
%these are for the bottom ports
z_shift = 2*3.2/2^11; %two cell inward
z_bot_str = -0.47919493 - z_shift;
z_bot_end = -0.48329875 - z_shift;

Nx = 25;
Nz = Nx;
Ny = 45;
xs = linspace(x_str,x_end,Nx);
zs_bot= linspace(z_bot_str,z_bot_end,Nz);
zs_top = -zs_bot;
%% start deploying the points
% start with top ports, from negative y to positive y
% then with bot ports,  from negative y to positive y
% 6 ports on each end
Pos_top_ports=zeros(Nx*Ny*6,3);
for i = 1:6
    %for each port
    %% 
    y_str = ys_str(i);
    ys = linspace(y_str,y_str+dy,Ny);
    for j = 1:Ny
        %do it for each y in ys
        for k = 1:Nx
            %for each point with x and z at a given y
            y = ys(j);
            x = xs(k);
            z = zs_top(k);
            ind = (i-1)*Nx*Ny+(j-1)*Nx + k;
            Pos_top_ports(ind,1) = x;
            Pos_top_ports(ind,2) = y;
            Pos_top_ports(ind,3) = z;
        end
    end
end
%%

figure;
scatter3(Pos_top_ports(:,1),Pos_top_ports(:,2),Pos_top_ports(:,3));
hold on;
xlabel("x")
ylabel("y")
zlabel("z")
trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), scaledVertices(:, 3),'FaceColor','r','FaceAlpha',0.3)
daspect([1 1 1])
%% 
Pos_bot_ports=zeros(Nx*Ny*6,3);
for i = 1:6
    %for each port
    %% 
    y_str = ys_str(i);
    ys = linspace(y_str,y_str+dy,Ny);
    for j = 1:Ny
        %do it for each y in ys
        for k = 1:Nx
            %for each point with x and z at a given y
            y = ys(j);
            x = xs(k);
            z = zs_bot(k);
            ind = (i-1)*Nx*Ny+(j-1)*Nx + k;
            Pos_bot_ports(ind,1) = x;
            Pos_bot_ports(ind,2) = y;
            Pos_bot_ports(ind,3) = z;
        end
    end
end

%% plot to check port data
Pos_ports=[Pos_top_ports;Pos_bot_ports];
figure;
scatter3(Pos_ports(:,1),Pos_ports(:,2),Pos_ports(:,3)); hold on;
xlabel("x")
ylabel("y")
zlabel("z")
trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), scaledVertices(:, 3),'FaceColor','r','FaceAlpha',0.3)
daspect([1 1 1])
%% output
writematrix(Pos_ports,'ports_probes_xyz_in_no_header.txt',Delimiter=" ")

%% create additional stack ports input (on equal z plane)

%% new ports output
%% use the same bounds but put it more on the inseide!
z_top_str = 0.48329875 + z_shift;
z_top_end = 0.53133694 - z_shift;

Npos = 5;

z_top = linspace(z_top_str,z_top_end,Npos);
z_bot = fliplr(linspace(-z_top_str,-z_top_end,Npos));
z_pos = [z_bot,z_top]

Pos_ports_new = zeros(Nx*Ny*6*2*Npos,3);
inds = zeros(Nx*Ny*6*2*Npos,1);
for i = 1:6
    y_str = ys_str(i);
    ys = linspace(y_str,y_str+dy,Ny);
    for Nlayer = 1:2*Npos
        for indy = 1:Ny
            for indx = 1:Nx
               %this is for each pair of ports go though
               %determine x,y,z coords
               x = xs(indx);
               y = ys(indy);
               z = z_pos(Nlayer);

               ind = (i-1)*Nx*Ny*2*Npos + (Nlayer-1)*(Nx*Ny) + (indy-1)*Nx + indx;
               inds(ind)=1;
               Pos_ports_new(ind,1) = x;
               Pos_ports_new(ind,2) = y;
               Pos_ports_new(ind,3) = z;
            end
        end
    end
end

%% validate and save new

figure;
scatter3(Pos_ports_new(:,1),Pos_ports_new(:,2),Pos_ports_new(:,3)); hold on;
xlabel("x")
ylabel("y")
zlabel("z")
trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), scaledVertices(:, 3),'FaceColor','r','FaceAlpha',0.3)
daspect([1 1 1])

writematrix(Pos_ports_new,'ports_probes_xyz_new_no_header.txt',Delimiter=" ")
%% TODO: probe locations at outlet 
%for mass flux estimation data and thrust
delta=3.2/256;
%total number of points: 6279 161x39
ys = -1+delta:delta:1-delta;
Ny = length(ys)
zs = 0:delta:0.49814422;
zs = [-flip(zs(2:end)),zs];
Nz = length(zs)
Pos_outlet = zeros(Ny*Nz,3);
for i = 1:Ny
    for j = 1:Nz
        %iterate through z then y
        ind = (i-1)*Nz + j;
        y =ys(i);
        z = zs(j);
        Pos_outlet(ind,1) = -delta; %x = -0.025/h one cell inwards
        Pos_outlet(ind,2) = y;
        Pos_outlet(ind,3) = z;
    end 
end
%% 
figure;
scatter3(Pos_outlet(:,1),Pos_outlet(:,2),Pos_outlet(:,3)), hold on
trisurf(fv.ConnectivityList, scaledVertices(:, 1), scaledVertices(:, 2), scaledVertices(:, 3),'FaceColor','r','FaceAlpha',0.3)
daspect([1 1 1])
%% output 
writematrix(Pos_outlet,"outlet_probes_xyz_in_no_header.txt",Delimiter=" ")
%% set centerline data
delta = 0.05;
x_c = -2:delta:30;
Nc = length(x_c);
Pos_center=zeros(Nc,3);
for i = 1:Nc
    Pos_center(i,1)=x_c(i);
    Pos_center(i,2)=0; %y-coord = 0 at centerline
    Pos_center(i,3)=0; %z-coord = 0 at cetnerline
end
writematrix(Pos_center,"center_probes_xyz_no_header.txt",Delimiter=" ")

%% set up axial line probes at multiple locations
clear;
delta = 1;
x_c = 0:delta:80;
% setting up nine lines (one in center, 4 at middle of the lips and 4 at
% the corners)
% organized as below (looking towards upstream) 
% 1 | 2 | 3
% 4 | 5 | 6
% 7 | 8 | 9
y_locs = [-1,0,1,-1,0,1,-1,0,1];
z_locs = [0.5,0.5,0.5,0,0,0,-0.5,-0.5,-0.5];
Nlines = length(y_locs); %number of lines
Np     = length(x_c); %probes per line
assert(Nlines==length(z_locs));
Pos_axials = zeros(Np*Nlines,3);
for i = 1:Nlines
    y = y_locs(i);
    z = z_locs(i);
    for j = 1:Np
        %iterate thru each line
        ind = (i-1)*Np+j;
        Pos_axials(ind,1) = x_c(j); % x-pos
        Pos_axials(ind,2) = y; %y-pos
        Pos_axials(ind,3) = z; %z-pos
    end
end
figure;
scatter3(Pos_axials(:,1),Pos_axials(:,2),Pos_axials(:,3))
writematrix(Pos_axials,"Axial_probes_xyz_no_header.txt",Delimiter=" ")
%% set a shearlayer line analysis (minor plane)
delta = 0.05;
x_str = 0;
x_end = 20;
z_str = 0.5;
z_end = 1.2;
xz_slope = (z_end-z_str)/(x_end-x_str);
theta = atan(xz_slope)
dx = delta * cos(theta)
dz = delta * sin(theta)
x = x_str:dx:x_end;
z = z_str:dz:z_end;
assert(length(x) == length(z))
Ns = length(x)
Pos_shear=zeros(Ns,3);
for i = 1:Ns
    Pos_shear(i,1) = x(i); %x-coord
    Pos_shear(i,2) = 0;    %y-coord
    Pos_shear(i,3) = z(i);  %z-coord
end
%% set up a signle port line analysis
delta = 3.2/256; 
L = 0.6;
zp_bot = -L:delta:0;
zp_top = flip(-zp_bot);
Np = length(zp_bot);
ys = [0.1428632,0.42856842,0.71427368];
ys = [flip(-ys),ys];
Ny = length(ys);
Pos_port_top = zeros(Np*Ny,3);
Pos_port_bot = zeros(Np*Ny,3);

for j = 1:length(ys)
    for i = 1:Np
        ind = (j-1)*Np + i;
        Pos_port_top(ind,1) = -0.2524895; %x-coord
        Pos_port_top(ind,2) = ys(j); %y-coord
        Pos_port_top(ind,3) = zp_top(i); %z-coord

        Pos_port_bot(ind,1) = -0.2524895; %x-coord
        Pos_port_bot(ind,2) = ys(j); %y-coord
        Pos_port_bot(ind,3) = zp_bot(i); %z-coord

    end
end

Pos_ports = [Pos_port_bot; Pos_port_top];
%%
figure; scatter3(Pos_ports(:,1),Pos_ports(:,2),Pos_ports(:,3));
%% 
writematrix(Pos_ports,"port_line_probes_xyz_no_header.txt",Delimiter=" ")



