%% main_L1_GRID.m
% Synthesize distributed controllers (L1, Oracle, Spatial Regret)
% for the 16-node grid system (4x4) used in the paper.
%
% Workflow:
%   1) Build the 4x4 grid topology used in the paper
%   2) Synthesize Oracle controller (SLS, L1 problem)
%   3) Synthesize Spatial Regret controller (L1, SLS)
%   4) Synthesize baseline L1 controller (SLS)
%   5) Compare L1 norms and plot results
%
% See also: main_SDP_GRID.m, plots_for_L1.m

clear all; close all; clc;
addpath('./Functions_SpRegret');  % Add path to the folder with auxiliary functions
rng(100);  % Set random seed for reproducibility

%% Parameters for the synthesis optimization
N_tf = 40;     N_FIR_closed_loop = 40;
NUMBER_OF_SAMPLING_POINTS = 1;
method_oracles = 'sls';
problem_oracles = 'l1';

reduce_the_model_after = false;

%% Defining Graph structure and Oracle structure
shape_of_the_mesh = [4,4];
n_agents = prod(shape_of_the_mesh);

% Use the exact 16-node grid adjacency used in the paper
% (Uncomment the next line to try a random grid instead)
% Adjacency_matrix = random_mesh_grid(shape_of_the_mesh, 0.25);
Adjacency_matrix = generate_grid_16_adjacency();

Graph = digraph(Adjacency_matrix~=0);
delays_matrix = distances(Graph);

% fprintf("Network topology generated: %d agents in a %dx%d grid\n", ...
%     n_agents, shape_of_the_mesh(1), shape_of_the_mesh(2));

%% Defining and creating the folders where to save results
savingfolder_results = './results/'+string(n_agents)+'_agents/';
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

%% Defining the oracle delays' structure

% Everybody to 4
oracle_graph = Adjacency_matrix;
oracle_graph([2,3,7] , 4) = 1;
oracle_delays = distances(digraph(oracle_graph~=0));




%% Oracle
fprintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ORACLE Synthesis          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");

N_tf_oracle = N_tf;
options_oracle = get_default_options('N_tf', N_tf,...
    'method', 'sls', ...
    'verbose', 0, ...
    'solver','gurobi', ... 
    'N_FIR_closed_loop', N_FIR_closed_loop);

disp(   ">Oracle N_tf    : " +num2str(options_oracle.N_tf))
disp(   ">Oracle Method  : " + options_oracle.method );
disp(   ">Oracle Problem : " + problem_oracles);


tStart = tic;
[K_oracle, Phis_oracle, objective_oracle] = calculus_distributed(sys,problem_oracles, oracle_delays, options_oracle);
% disp("Optimization result: "+ num2str(H2_norm_oracle));
tEnd = toc(tStart);
fprintf("END Oracle     (time:%.4f s)\n",tEnd);

% Computing the LFT tf
lft_oracle = my_lft(sys,K_oracle);
if reduce_the_model_after
    disp("Reducing the order of lft_oracle...") %#ok<UNRCH>
    tStart = tic;
    lft_oracle = minreal(lft_oracle);
    tEnd = toc(tStart);
    fprintf("Finished the order minimization (time:%.4f s)\n",tEnd);
end


%% Spatial Regret with SAMPLING
fprintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Spatial Regret Synthesis         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
options_spatial_regret = get_default_options('N_tf', N_tf,...
    'method', 'sls', ...
    'verbose', 0, ...
    'solver','gurobi', ... 
    'N_FIR_closed_loop', N_FIR_closed_loop);


tStart = tic;
[K_sp_reg_Sampled, Phis_spreg, spreg_l1] = calculus_spatial_regret_L1(sys, Phis_oracle, delays_matrix, options_spatial_regret);
tEnd = toc(tStart);
fprintf("END Spatial Regret (time:%.4f s)\n",tEnd);
% fprintf("------------------------------------------------------\n");    
lft_spreg = my_lft(sys,K_sp_reg_Sampled);
if reduce_the_model_after
    disp("Reducing the order of lft_spreg...") %#ok<UNRCH>
    tStart = tic;
    lft_spreg = minreal(lft_spreg);
    tEnd = toc(tStart);
    fprintf("Finished the order minimization (time:%.4f s)\n",tEnd);
end


%% SLS L1
fprintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ L1 Distributed Synthesis         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
options_l1_sls = get_default_options('N_tf', N_tf,...
    'method', 'sls', ...
    'verbose', 0, ...
    'solver','gurobi', ... 
    'N_FIR_closed_loop', N_FIR_closed_loop);


tStart = tic;
[K_l1, Phis_l1, L1_norm] = calculus_distributed(sys,'l1',delays_matrix, options_l1_sls);
tEnd = toc(tStart);
fprintf("END L1 SLS  (time:%.4f s)\n",tEnd);
fprintf("------------------------------------------------------\n");

lft_l1 = my_lft(sys,K_l1);

if reduce_the_model_after
    fprintf("Reducing the order of lft_L1...\n"); %#ok<UNRCH>
    tStart = tic;
    lft_l1 = minreal(lft_l1);
    tEnd = toc(tStart);
    fprintf("Finished the order minimization (time:%.4f s)\n",tEnd);
end

fprintf("L1 norm: %.4f\n", L1_norm);
fprintf("------------------------------------------------------\n");


%% Saving
save(savingfolder_results + "results_with_NFIR_"+num2str(N_tf)+"_N_Cl_" + num2str(N_FIR_closed_loop)+ ".mat");

%% Comparing the results performance
fprintf("Comparing the results of the controllers...\n");

% Prepare storage for norms and labels
controllerNames = {};
L1_norms = [];
H2_norms = [];
Hinf_norms = [];

% Row 1: L1 controller
controllerNames{end+1} = 'L1 controller';
L1_norms(end+1)        = L1_norm;  % objective from L1 synthesis

% Now loop over all oracle and spatial regret controllers.
    % For the Oracle controller corresponding to iteration i:
    fprintf("    Oracle -->")
    L1_norm_oracle = L1norm_FIR(sys, Phis_oracle.mat ,max(options_oracle.N_FIR_closed_loop, N_tf));
    fprintf(" Done\n");
    controllerNames{end+1} = 'Oracle';
    L1_norms(end+1)        = L1_norm_oracle;
    
    fprintf("    SpReg -->")
    L1_norm_spreg = L1norm_FIR(sys, Phis_spreg.mat ,max(options_oracle.N_FIR_closed_loop, N_tf));
    fprintf(" Done\n");
    controllerNames{end+1} = 'SpRegret';
    L1_norms(end+1)        = L1_norm_spreg;


% Create and display table
disp("------------------------------------------------------");
disp("Norms table");
disp("------------------------------------------------------");
norm_table = table(L1_norms', ...
    'VariableNames', {'L1_Norm'}, 'RowNames', controllerNames);
disp(norm_table);

%% Plotting
save_plots_as_figs = true; % Save plots also as .fig files (for interactive editing)
plots_for_L1;

close all;

chosen_channel = 4;
Sim_one_specific_disturbance_L1;


%% Computing the solution in a distributed fashion

options_distribut = get_default_options('N_tf', N_tf,...
    'method', 'sls', ...
    'verbose', 0, ...
    'solver','gurobi', ... 
    'N_FIR_closed_loop', N_tf);

    % Phi_warm_start_mat = struct;
    % Phi_warm_start_mat.R = Rmat;
    % Phi_warm_start_mat.N = Nmat;
    % Phi_warm_start_mat.M = Mmat;
    % Phi_warm_start_mat.L = Lmat;
    % options_distribut.warmstart =  true;
    % options_distribut.warmStartSolution = Phi_warm_start_mat;
    % options_distribut.warmstart =  true;
    options_distribut.warmStartSolution = Phis_l1.mat;

    options_distribut.distributed_optimization = true;

    options_distribut.rho = 0.01;      % Initial penalty parameter 
    options_distribut.max_ADMM_iter = 1e4;
    options_distribut.tol_ADMM = 1e-3; % Tolerance for ADMM convergence
    options_distribut.tol_gr = 1e-3; % Tolerance for convergence
    options_distribut.max_iter_gr = 100;
    options_distribut.gamma_upper = 250; % Upper bound for the bisection method
    options_distribut.gamma_lower = 0.01; % Lower bound for the bisection method

    % Set the row and column indices for the decision variables
    options_distribut.N_row_indices = 16;
    options_distribut.N_col_indices = 1;

    [K_opt_distrib, Phis_spreg_distrib, objective_distrib]  = calculus_spatial_regret_L1(sys, Phis_oracle, delays_matrix, options_distribut);