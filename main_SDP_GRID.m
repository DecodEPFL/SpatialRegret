%% main_SDP_GRID.m
% Synthesize distributed controllers (H2, Hinf, Oracle, Spatial Regret)
% for a 5-node linear system using Youla parameterization and SDP methods.
% This script has been used to generate the results in Section IV.A of the paper.
%
% Workflow:
%   1) Build a 5-agent linear chain topology
%   2) Synthesize Oracle controller (Hinf, full information sharing to last agent)
%   3) Synthesize Spatial Regret controller (minimizes worst-case regret w.r.t. oracle)
%   4) Synthesize baseline H2 and Hinf controllers
%   5) Compare performance using H2, Hinf, and Spatial Regret norms
%   6) Simulate impulse responses and generate plots
%
% See also: main_L1_GRID.m, plots_for_SDP.m

clear all; close all; clc;
addpath('./Functions_SpRegret');  % Add path to the folder with auxiliary functions
rng(100);  % Set random seed for reproducibility

%% Parameters for the synthesis optimization
N_tf = 20;
NUMBER_OF_SAMPLING_POINTS = 4000;

%% Defining Graph structure and Oracle structure:
shape_of_the_mesh = [1,5];

n_agents = prod(shape_of_the_mesh);

Adjacency_matrix = eye(n_agents);
Adjacency_matrix = Adjacency_matrix + diag(ones(n_agents-1,1),1) + diag(ones(n_agents-1,1),-1);


Graph = digraph(Adjacency_matrix~=0);
delays_matrix = distances(Graph);

%% Defining and creating the folders where to save results
savingfolder_results = './results/' + string(n_agents) + '_agents/';
if ~exist(savingfolder_results, 'dir')
    mkdir(savingfolder_results)
end
savingfolder_figures = './figures/'+string(n_agents)+'_agents/';
if ~exist(savingfolder_figures, 'dir')
    mkdir(savingfolder_figures)
end



%% Defining the plant model
Ts = 0.1;  % Sampling period
disp("Generating the plant...");
is_state_feedback = false;
plant = generate_plant_homogeneous(Adjacency_matrix, Ts,is_state_feedback);

sys.A = plant.A;
sys.B2 = plant.B;
sys.C2 = plant.C;
sys.D22 = plant.D;
sys.plant = plant;



%% Defining weight Matrices
sys.n = size(sys.A, 1);    % Order of the system: state dimension
sys.m = size(sys.B2, 2);   % Number of input channels
sys.p = size(sys.C2,1);    % Number of output channels

Q_weight = 1*eye(sys.n);  %Weights for the states
R_weight = 1*eye(sys.m);  %Weights for the inputs

% Performance matrices such that z_t' * z_t := x_t' * Q * x_t  +  u_t' * R * u_t
sys.C1 = [sqrtm(Q_weight);
    zeros(sys.m, sys.n)];

sys.D12= [zeros(sys.n,sys.m);
    sqrtm(R_weight)];

sys.n_z = size(sys.C1,1); % Number of performance outputs

%% Defining disturbance Matrices

sys.n_w = n_agents; % Number of disturbances
sys.B1 = kron(eye(n_agents),[0;1]); % Disturbances hitting the second state of each agent
sys.D11 = zeros(sys.n_z, sys.n_w);
sys.D21 = eye(sys.p,sys.n_w);




%% Defining Doubly-Coprime Factorization
[F, L] = calculus_F_and_L(sys, Adjacency_matrix); % Computes the sparse F and L matrices
sys.F = F;
sys.L = L;

[P11,P12,P21] = coprime_factorization(sys); % Computes the P11, P12 and P21 matrices
sys.P11 = P11;
sys.P12 = P12;
sys.P21 = P21;

%% Defining the oracle supergraph

% Everbody to last
oracle_graph = Adjacency_matrix;
oracle_graph(:, end) = 1;
oracle_delays = distances(digraph(oracle_graph~=0));


%% Oracle
fprintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ORACLE SYNTHESIS          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");

N_tf_oracle = N_tf;
options_oracle =  get_default_options('N_tf', N_tf_oracle,...
    'method','youla');

disp(   ">Oracle N_tf    : " +num2str(options_oracle.N_tf))
disp(   ">Oracle Method  : " + options_oracle.method );


tStart = tic;
[K_h2_oracle, Q_oracle, H2_norm_oracle] = calculus_distributed(sys, 'hinf', oracle_delays, options_oracle);
% disp("Optimization result: "+ num2str(H2_norm_oracle));
tEnd = toc(tStart);
fprintf("END Oracle.  Opt val: %.4f   (time:%.4f s)\n",H2_norm_oracle, tEnd);

% Computing the Oracle LFT
lft_oracle = minreal(sys.P11 - sys.P12*Q_oracle*sys.P21);

%% Spatial Regret with SAMPLING
fprintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Spatial Regret Synthesis         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
options_spatial_regret_sampled = get_default_options('number_points', NUMBER_OF_SAMPLING_POINTS, 'N_tf', N_tf,...
    'method', 'tf_sampled');

tStart = tic;
[K_sp_reg_Sampled, Q_spreg, final_spreg_sampled] = calculus_spatial_regret(sys, lft_oracle, delays_matrix, options_spatial_regret_sampled);
tEnd = toc(tStart);
fprintf("END Spatial Regret.  Opt val: %.4f   (time:%.4f s)\n",final_spreg_sampled, tEnd);
fprintf("------------------------------------------------------\n");

lft_spreg = minreal(sys.P11 - sys.P12*Q_spreg*sys.P21);


%% H2 Distributed Youla
fprintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ H2 Distributed Synthesis         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");

options_h2_youla = get_default_options('number_points', NUMBER_OF_SAMPLING_POINTS, 'N_tf', N_tf,...
    'method', 'youla');

tStart = tic;
[K_h2_dist_youla, Q_H2, H2_norm_youla] = calculus_distributed(sys,'h2',delays_matrix, options_h2_youla);
tEnd = toc(tStart);
fprintf("END Youla H2.  Opt val: %.4f   (time:%.4f s)\n",H2_norm_youla, tEnd);
fprintf("------------------------------------------------------\n");

lft_h2 = minreal(sys.P11 - sys.P12*Q_H2*sys.P21);


%% HINF
fprintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HINF Distributed Synthesis         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");

options_hinf_youla = get_default_options('number_points', NUMBER_OF_SAMPLING_POINTS, 'N_tf', N_tf,...
    'method', 'youla');

tStart = tic;
[~, Q_Hinf, ~] = calculus_distributed(sys,'hinf',delays_matrix, options_hinf_youla);
tEnd = toc(tStart);
fprintf("END Youla Hinf  (time:%.4f s)\n\n",tEnd);

lft_hinf = minreal(sys.P11 - sys.P12*Q_Hinf*sys.P21);



%% H2 Centralized Controllers (for baseline comparison)
fprintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ H2 Centralized synthesis         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
tStart = tic;
total_sys = ss(sys.A,[sys.B1, sys.B2], [sys.C1; sys.C2],[sys.D11, sys.D12;sys.D21, sys.D22],Ts);
[K_h2_mat,CL_h2_mat,h2_norm_mat] = h2syn(total_sys,sys.p,sys.m);
tEnd = toc(tStart);
fprintf("END H2 Centralized.  Opt val: %.4f   (time:%.4f s)\n",h2_norm_mat, tEnd);
fprintf("------------------------------------------------------\n");

%% Saving
save(savingfolder_results + "results_SDP_with_NFIR_"+num2str(N_tf)+".mat");

%% Plotting
save_plots_as_figs = true; % Set to true to save plots as .fig files

plots_for_SDP;

chosen_channel = 5; % disturbance channel to affect (0 for all);
Sim_one_specific_disturbance_SDP;


%% Final Table with norms
% Compute the three norms for each controller
disp("Final Table with norms");
disp("------------------------------------------------------");
h2_vals = [ norm(lft_h2,2)^2;
    norm(lft_spreg,2)^2;
    norm(lft_hinf,2)^2];

hinf_vals = [ hinfnorm(lft_h2)^2;
    hinfnorm(lft_spreg)^2;
    hinfnorm(lft_hinf)^2
    ];

spreg_vals = [ spregnorm(lft_h2, lft_oracle, 3000);
    spregnorm(lft_spreg, lft_oracle, 3000);
    spregnorm(lft_hinf, lft_oracle, 3000)];

% Build a MATLAB table with row names
controller_names = { 'lft_h2'
    'lft_spreg',
    'lft_hinf' };

T = table(h2_vals, hinf_vals, spreg_vals, ...
    'RowNames', controller_names, ...
    'VariableNames', {'H2^2', 'Hinf^2', 'SpRegNorm'});

disp(T)