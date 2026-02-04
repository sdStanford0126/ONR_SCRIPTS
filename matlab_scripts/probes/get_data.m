% get data and store in matrix
% this outputs u-velocity and temperature
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE THIS BASED ON YOUR SAMPLING MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy = 0.05;
dz = 0.05;
x_samp = [0,0.5,1.8,3.5,6.5,13,19,24,31,40];
y_samp = [-10:dy:10];
z_samp = [-10:dz:10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PROBE FILE INFO 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prob_prefix = ''



filename1 = strcat(['probe_extended.00820002.pcd']);
filename2 = strcat(['probe_extended.pxyz']);
fprintf(filename1)
fprintf(filename2)
data = importdata(filename1);
xyz = importdata(filename2);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE THIS BASED ON VARIABLES YOU SAMPLED AND HOW THEY ARE STORED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indices = xyz.data(:,1)+1;
x = xyz.data(:,2);
y = xyz.data(:,3);
z = xyz.data(:,4);
u = data.data(:,2);
T = data.data(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_mat = zeros(length(y_samp), length(x_samp), length(z_samp));
T_mat = zeros(length(y_samp), length(x_samp), length(z_samp));
u_mat_x = zeros(length(y_samp), length(x_samp));
T_mat_x = zeros(length(y_samp), length(x_samp));

fprintf('Done setting things up.. Enter looping')

parpool('local', 100);   
for i = 1:length(z_samp)
    map = find(abs(z - z_samp(i)) < 0.001);
    indices_new = indices(map);
    u_new = u(map);
    T_new = T(map);

    parfor j = 1:length(x_samp)
        u_tmp = zeros(length(y_samp),1);
        T_tmp = zeros(length(y_samp),1);

	    for k = 1:length(y_samp)
    
		    counter = (i-1)*length(x_samp)*length(y_samp) + (j-1)*length(y_samp) + k;
		    idx = find(abs(indices_new - counter) < 0.1);
		    u_tmp(k) = u_new(idx);
		    T_tmp(k) = T_new(idx);
	    end
	    u_mat_x(:,j) = u_tmp;
	    T_mat_x(:,j) = T_tmp;
    end

    u_mat(:,:,i) = u_mat_x;
    T_mat(:,:,i) = T_mat_x;

end

fprintf('SAVING DATA\n')
filename3 = strcat(['data.mat']);
save(filename3, 'u_mat', 'T_mat', 'x_samp', 'y_samp', 'z_samp', '-v7.3');

delete(gcp('nocreate'));
