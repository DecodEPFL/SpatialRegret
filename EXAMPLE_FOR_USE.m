%% Example 2: Spatial Regret Comparison
% This example demonstrates how to compute and compare spatial regret
% for different controller architectures.

clear all; close all; clc;

% Add path to functions
addpath('./Functions_SpRegret');

%% Step 1: Setup (same as Example 1)
n_agents = 4;
Adjacency_matrix = eye(n_agents) + ...
                   diag(ones(n_agents-1,1), 1) + ...
                   diag(ones(n_agents-1,1), -1);

Graph = digraph(Adjacency_matrix ~= 0);
delays_matrix = distances(Graph);

Ts = 0.1;
plant = generate_plant_homogeneous(Adjacency_matrix, Ts, false);

sys.A = plant.A;
sys.B2 = plant.B;
sys.C2 = plant.C;
sys.D22 = plant.D;
sys.plant = plant;

sys.n = size(sys.A, 1);
sys.m = size(sys.B2, 2);
sys.p = size(sys.C2, 1);

Q_weight = eye(sys.n);
R_weight = eye(sys.m);

sys.C1 = [sqrtm(Q_weight); zeros(sys.m, sys.n)];
sys.D12 = [zeros(sys.n, sys.m); sqrtm(R_weight)];
sys.n_z = size(sys.C1, 1);

sys.n_w = n_agents;
sys.B1 = kron(eye(n_agents), [0; 1]);
sys.D11 = zeros(sys.n_z, sys.n_w);
sys.D21 = eye(sys.p, sys.n_w);

[F, L] = calculus_F_and_L(sys, Adjacency_matrix);
sys.F = F;
sys.L = L;

[P11, P12, P21] = coprime_factorization(sys);
sys.P11 = P11;
sys.P12 = P12;
sys.P21 = P21;

%% Step 2: Define Oracle Structure
% Oracle: All agents can share information with the last agent
fprintf('--- Defining Oracle Information Structure ---\n');

oracle_delays = delays_matrix;
oracle_delays(:, end) = 0;  % All agents can instantly communicate with last agent

fprintf('Standard communication delays:\n');
disp(delays_matrix);
fprintf('\nOracle communication delays (enhanced info sharing):\n');
disp(oracle_delays);

%% Step 3: Synthesis Options
options_oracle = get_default_options(...
    'number_points', 2000, ...
    'N_tf', 20, ...
    'method', 'youla');

options_spreg = get_default_options(...
    'number_points', 2000, ...
    'N_tf', 20, ...
    'method', 'tf_sampled');

%% Step 4: Synthesize Oracle Controller
fprintf('\n--- Synthesizing Oracle Controller ---\n');
tic;
[K_oracle, Q_oracle, H2_oracle] = calculus_distributed(sys, 'hinf', ...
                                                       oracle_delays, ...
                                                       options_oracle);
oracle_time = toc;

lft_oracle = sys.P11 - sys.P12 * Q_oracle * sys.P21;
lft_oracle = minreal(lft_oracle);

fprintf('Oracle synthesis completed in %.2f seconds\n', oracle_time);
fprintf('Oracle H-infinity norm: %.4f\n', H2_oracle);

%% Step 5: Synthesize Spatial Regret Controller
fprintf('\n--- Synthesizing Spatial Regret Controller ---\n');
tic;
[K_spreg, Q_spreg, spreg_cost] = calculus_spatial_regret(sys, ...
                                                         lft_oracle, ...
                                                         delays_matrix, ...
                                                         options_spreg);
spreg_time = toc;

lft_spreg = sys.P11 - sys.P12 * Q_spreg * sys.P21;
lft_spreg = minreal(lft_spreg);

fprintf('Spatial regret synthesis completed in %.2f seconds\n', spreg_time);
fprintf('Spatial regret cost: %.4f\n', spreg_cost);

%% Step 6: Synthesize Standard H2 Controller (for comparison)
fprintf('\n--- Synthesizing Standard H2 Controller ---\n');
tic;
[K_h2, Q_h2, H2_h2] = calculus_distributed(sys, 'h2', ...
                                          delays_matrix, ...
                                          options_oracle);
h2_time = toc;

lft_h2 = sys.P11 - sys.P12 * Q_h2 * sys.P21;
lft_h2 = minreal(lft_h2);

fprintf('H2 synthesis completed in %.2f seconds\n', h2_time);
fprintf('H2 norm: %.4f\n', H2_h2);

%% Step 7: Compute Performance Metrics
fprintf('\n--- Computing Performance Metrics ---\n');

% H2 norms
h2_oracle = norm(lft_oracle, 2);
h2_spreg = norm(lft_spreg, 2);
h2_h2ctrl = norm(lft_h2, 2);

% H-infinity norms
hinf_oracle = hinfnorm(lft_oracle);
hinf_spreg = hinfnorm(lft_spreg);
hinf_h2ctrl = hinfnorm(lft_h2);

% Spatial regret norms (requires custom function)
    N_points = 2000;
    spregnorm_oracle = 0;  % Oracle compared to itself
    spregnorm_spreg = spregnorm(lft_spreg, lft_oracle, N_points);
    spregnorm_h2ctrl = spregnorm(lft_h2, lft_oracle, N_points);


%% Step 8: Display Results Table
fprintf('\n=== PERFORMANCE COMPARISON ===\n\n');

controller_names = {'Oracle (Hinf)', 'Spatial Regret', 'Standard H2'};
h2_norms = [h2_oracle; h2_spreg; h2_h2ctrl];
hinf_norms = [hinf_oracle; hinf_spreg; hinf_h2ctrl];
spreg_norms = [spregnorm_oracle; spregnorm_spreg; spregnorm_h2ctrl];

results_table = table(h2_norms, hinf_norms, spreg_norms, ...
    'RowNames', controller_names, ...
    'VariableNames', {'H2_Norm', 'Hinf_Norm', 'SpRegret_Norm'});

disp(results_table);

%% Step 9: Frequency Response Comparison
fprintf('\n--- Generating Frequency Response Plots ---\n');

omega = logspace(-2, 2, 200);  % Frequency range

% Compute frequency responses
[mag_oracle, ~] = bode(lft_oracle, omega);
[mag_spreg, ~] = bode(lft_spreg, omega);
[mag_h2, ~] = bode(lft_h2, omega);

% Convert to arrays and compute max singular values
mag_oracle = squeeze(mag_oracle);
mag_spreg = squeeze(mag_spreg);
mag_h2 = squeeze(mag_h2);

% Maximum singular value at each frequency
sv_oracle = zeros(1, length(omega));
sv_spreg = zeros(1, length(omega));
sv_h2 = zeros(1, length(omega));

for i = 1:length(omega)
    sv_oracle(i) = max(svd(mag_oracle(:,:,i)));
    sv_spreg(i) = max(svd(mag_spreg(:,:,i)));
    sv_h2(i) = max(svd(mag_h2(:,:,i)));
end

% Plot
figure('Name', 'Frequency Response Comparison');
semilogx(omega, 20*log10(sv_oracle), 'b-', 'LineWidth', 2, 'DisplayName', 'Oracle');
hold on;
semilogx(omega, 20*log10(sv_spreg), 'r--', 'LineWidth', 2, 'DisplayName', 'Spatial Regret');
semilogx(omega, 20*log10(sv_h2), 'g-.', 'LineWidth', 2, 'DisplayName', 'Standard H2');
hold off;
grid on;
xlabel('Frequency (rad/s)');
ylabel('Max Singular Value (dB)');
title('Closed-Loop Frequency Response');
legend('Location', 'best');

%% Summary
fprintf('\n=== KEY INSIGHTS ===\n');
fprintf('1. Oracle controller has access to enhanced information structure\n');
fprintf('2. Spatial Regret controller optimizes worst-case regret w.r.t. oracle\n');
fprintf('3. Standard H2 only optimizes average performance\n\n');

fprintf('Performance degradation (vs Oracle H2 norm):\n');
fprintf('  Spatial Regret: %.2f%%\n', 100*(h2_spreg - h2_oracle)/h2_oracle);
fprintf('  Standard H2:    %.2f%%\n\n', 100*(h2_h2ctrl - h2_oracle)/h2_oracle);

fprintf('Example completed successfully!\n');
