% SPOD ANALYSIS
% Modified by Steven Dai baded on script by Olivia Martin
% This SPOD analysis computes SPOD using variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_dir = "./baseline_new/130M/";
data_name = "pressure_xz_plane.mat";
load(input_dir + data_name)

p = P(:,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = 1.4;
numT = size(p,1); % Number of time samples
num_vars = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE WEIGHT MATRIX DIMENSIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numSnaps = numT/2+1;
numBlocks = numT*3/numSnaps - 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE SPOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weight = [];
window = hann(numSnaps);
noverlap = []; %numSnaps/2 - 2;
dt = 0.15;
opts.savefft    = true;             % save FFT blocks insteasd of keeping them in memory
opts.savedir    = 'spod_results_test/';
opts.mean       = 'blockwise';

[L,P,f,Lc,A] = spod(p, window,weight,noverlap, dt, opts);

save('spod_results_test/exp_coeff.mat', 'A');

