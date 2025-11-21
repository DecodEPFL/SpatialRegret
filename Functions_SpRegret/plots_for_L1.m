% clear all; close all; clc;
% addpath('./Functions_SpRegret');  % Add path to the folder with auxiliary functions
%% Plots' parameters

% Define the different parameter for the plots
line_width = 1.5;
saving_fig_file = true;
alpha_L1 = 1.0;
alpha_Spreg_L1 = 0.5;
plot_params = struct;
plot_params.half_frequencies = true; % Use half frequencies (0 to pi)
plot_params.plot_the_oracle = true;

plot_params.N_FREQ_POINTS = 1500;
plot_params.line_width_l1 = 2;
plot_params.line_width_spreg = 2;
plot_params.line_width_oracle = 1*2;
plot_params.with_logy = false;
plot_params.Position = [100,100,556,272];

plot_params.linestyle_l1 = '-';
plot_params.linestyle_spreg = '-';
plot_params.linestyle_oracle = '--';

plot_params.color_l1    = [0.2157, 0.4941, 0.7225];   % blue
plot_params.color_spreg = [0.8941, 0.1019, 0.1098];   % red
% plot_params.color_oracle = [0.8500 0.3250 0.0980];   % orange;
plot_params.color_oracle = 'k';   % orange;


if ~exist(savingfolder_figures, 'dir')
    mkdir(savingfolder_figures)
end

%% Frequency Response Plots

Phis_L1_mat    = compute_FIR_Matrix(sys, Phis_l1.mat, N_tf);
Phis_Spreg_mat = compute_FIR_Matrix(sys, Phis_spreg.mat, N_tf);
Phis_Oracle_mat = compute_FIR_Matrix(sys, Phis_oracle.mat, N_tf);


for input_channel = 1:n_agents
    plot_l1_channel_i = draw_L1_freq_plots_with_FIR(Phis_Spreg_mat, Phis_Oracle_mat, Phis_L1_mat, plot_params, input_channel);
    exportgraphics(plot_l1_channel_i, savingfolder_figures + "L1_freq_plot_on_channel_" + num2str(input_channel) + ".pdf", 'ContentType', 'vector');
    if saving_fig_file
        saveas(plot_l1_channel_i, savingfolder_figures + "L1_freq_plot_on_channel_" + num2str(input_channel) + ".fig");
    end
end


function [final_figure] = draw_L1_freq_plots(lft_spreg, lft_oracle, lft_l1, plot_params, input_channel)
if nargin < 5
    input_channel = 0;
end

N_FREQ_POINTS = plot_params.N_FREQ_POINTS;
frequency_vector = linspace(-pi, pi, N_FREQ_POINTS);

% Initialize cost arrays for the different controllers
cost_l1      = zeros(1, N_FREQ_POINTS);
cost_spreg   = zeros(1, N_FREQ_POINTS);
cost_oracle  = zeros(1, N_FREQ_POINTS);
disp("DEBUG: Starting to compute the costs");
if input_channel > 0
    lft_l1_jth_channel     = lft_l1(:,input_channel);
    lft_spreg_jth_channel  = lft_spreg(:,input_channel);
    lft_oracle_jth_channel = lft_oracle(:,input_channel);

    for k = 1:N_FREQ_POINTS
        z_k = exp(1i * frequency_vector(k));

        Tzw_l1     = evalfr(lft_l1_jth_channel, z_k);
        Tzw_spreg  = evalfr(lft_spreg_jth_channel, z_k);
        Tzw_oracle = evalfr(lft_oracle_jth_channel, z_k);

        cost_l1(k)     = norm(Tzw_l1, Inf);
        cost_spreg(k)  = norm(Tzw_spreg, Inf);
        cost_oracle(k) = norm(Tzw_oracle, Inf);

        % assert(cost_l1(k) == norm(Tzw_l1, Inf), "L1 cost is not equal to the infinity norm " + num2str(k));
    end
else
    for k = 1:N_FREQ_POINTS
        z_k = exp(1i * frequency_vector(k));

        Tzw_l1     = evalfr(lft_l1, z_k);
        Tzw_spreg  = evalfr(lft_spreg, z_k);
        Tzw_oracle = evalfr(lft_oracle, z_k);

        cost_l1(k)     = norm(Tzw_l1, Inf);
        cost_spreg(k)  = norm(Tzw_spreg, Inf);
        cost_oracle(k) = norm(Tzw_oracle, Inf);
    end
end

disp("DEBUG: Finished computing the costs");
color_l1 = plot_params.color_l1;
color_spreg = plot_params.color_spreg;
color_oracle = plot_params.color_oracle;

linestyle_l1 = plot_params.linestyle_l1;
linestyle_spreg = plot_params.linestyle_spreg;
linestyle_oracle = plot_params.linestyle_oracle;

line_width = plot_params.line_width;
with_logy = plot_params.with_logy;

final_figure = figure;
hold on;

if with_logy
    plot_l1 =  plot(frequency_vector, log10(cost_l1), 'Color', color_l1, 'LineStyle', linestyle_l1, 'LineWidth', line_width, 'DisplayName', "$\mathbf{K}^{\mathcal{L}_{1}}$");
    plot_oracle = plot(frequency_vector, log10(cost_oracle), 'Color', color_oracle, 'LineStyle', linestyle_oracle, 'LineWidth', line_width, 'DisplayName', "$\mathbf{\hat{K}}$");
    plot_spreg = plot(frequency_vector, log10(cost_spreg), 'Color', color_spreg, 'LineStyle', linestyle_spreg, 'LineWidth', line_width, 'DisplayName', "$\mathbf{K}^{{SR}}_{\mathcal{L}_{1}}$");
else
    plot_l1 =  plot(frequency_vector, cost_l1, 'Color', color_l1, 'LineStyle', linestyle_l1, 'LineWidth', line_width, 'DisplayName', "$\mathbf{K}^{\mathcal{L}_{1}}$");
    plot_oracle = plot(frequency_vector, cost_oracle, 'Color', color_oracle, 'LineStyle', linestyle_oracle, 'LineWidth', line_width, 'DisplayName', "$\mathbf{\hat{K}}$");
    plot_spreg =  plot(frequency_vector, cost_spreg, 'Color', color_spreg, 'LineStyle', linestyle_spreg, 'LineWidth', line_width, 'DisplayName', "$\mathbf{K}^{{SR}}_{\mathcal{L}_{1}}$");
end

legend([plot_spreg, plot_l1, plot_oracle], 'Interpreter','latex', 'FontSize',14);
xlabel('Frequency $\omega$', 'Interpreter','latex', 'FontSize',18);
if with_logy
    ylabel('$\log_{10}(||\mathbf{F}^{[:,' + string(input_channel) + ']}(e^{j\Omega})||)_{\infty}^2$', 'Interpreter','latex', 'FontSize',18);
else
    ylabel('$||\mathbf{F}^{[:,' + string(input_channel) + ']}(e^{j\Omega})||_{\infty}^2$', 'Interpreter','latex', 'FontSize',18);
end
grid minor;
xlim([-pi, pi]);
xticks([-pi, -pi/2, 0, pi/2, pi]);
xticklabels({'$-\pi$', '$-\frac{\pi}{2}$', '0', '$\frac{\pi}{2}$', '$\pi$'});
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
xaxisproperties.FontSize = 12;
if ~isempty(plot_params.Position)
    % Set the figure position if specified
    final_figure.Position = plot_params.Position;
end
hold off;
end


function [final_figure] = draw_L1_freq_plots_with_FIR(lft_spreg_matrix, lft_oracle_matrix, lft_l1_matrix, plot_params, input_channel, output_channel)
if nargin < 5
    input_channel = 0;
    output_channel = 0;
end

if nargin <6
    output_channel = 0;
end

half_frequencies = plot_params.half_frequencies;    
N_FREQ_POINTS = plot_params.N_FREQ_POINTS;
if half_frequencies
    frequency_vector = linspace(0, pi, N_FREQ_POINTS);
else
    frequency_vector = linspace(-pi, pi, N_FREQ_POINTS);
end

% Initialize cost arrays for the different controllers
cost_l1      = zeros(1, N_FREQ_POINTS);
cost_spreg   = zeros(1, N_FREQ_POINTS);
cost_oracle  = zeros(1, N_FREQ_POINTS);
% disp(" Input channel: " + num2str(input_channel));
if input_channel > 0 && output_channel == 0
    lft_l1_jth_channel     = lft_l1_matrix(:,input_channel,:);
    lft_spreg_jth_channel  = lft_spreg_matrix(:,input_channel,:);
    lft_oracle_jth_channel = lft_oracle_matrix(:,input_channel,:);

    for k = 1:N_FREQ_POINTS
        z_k = exp(1i * frequency_vector(k));

        Tzw_l1     = evalfr_FIR(lft_l1_jth_channel, z_k);
        Tzw_spreg  = evalfr_FIR(lft_spreg_jth_channel, z_k);
        Tzw_oracle = evalfr_FIR(lft_oracle_jth_channel, z_k);

        cost_l1(k)     = norm(Tzw_l1, Inf);
        cost_spreg(k)  = norm(Tzw_spreg, Inf);
        cost_oracle(k) = norm(Tzw_oracle, Inf);

        % assert(cost_l1(k) == norm(Tzw_l1, Inf), "L1 cost is not equal to the infinity norm " + num2str(k));
    end
elseif input_channel == 0 && output_channel == 0
    for k = 1:N_FREQ_POINTS
        z_k = exp(1i * frequency_vector(k));

        Tzw_l1     = evalfr_FIR(lft_l1_matrix, z_k);
        Tzw_spreg  = evalfr_FIR(lft_spreg_matrix, z_k);
        Tzw_oracle = evalfr_FIR(lft_oracle_matrix, z_k);

        cost_l1(k)     = norm(Tzw_l1, Inf);
        cost_spreg(k)  = norm(Tzw_spreg, Inf);
        cost_oracle(k) = norm(Tzw_oracle, Inf);
    end
elseif input_channel == 0 && output_channel > 0
    lft_l1_jth_channel     = lft_l1_matrix(output_channel,:,:);
    lft_spreg_jth_channel  = lft_spreg_matrix(output_channel,:,:);
    lft_oracle_jth_channel = lft_oracle_matrix(output_channel,:,:);

    for k = 1:N_FREQ_POINTS
        z_k = exp(1i * frequency_vector(k));

        Tzw_l1     = evalfr_FIR(lft_l1_jth_channel, z_k);
        Tzw_spreg  = evalfr_FIR(lft_spreg_jth_channel, z_k);
        Tzw_oracle = evalfr_FIR(lft_oracle_jth_channel, z_k);

        cost_l1(k)     = norm(Tzw_l1, Inf);
        cost_spreg(k)  = norm(Tzw_spreg, Inf);
        cost_oracle(k) = norm(Tzw_oracle, Inf);

        % assert(cost_l1(k) == norm(Tzw_l1, Inf), "L1 cost is not equal to the infinity norm " + num2str(k));
    end
else % input_channel > 0 && output_channel > 0
    lft_l1_jth_channel     = lft_l1_matrix(output_channel,input_channel,:);
    lft_spreg_jth_channel  = lft_spreg_matrix(output_channel,input_channel,:);
    lft_oracle_jth_channel = lft_oracle_matrix(output_channel,input_channel,:);

    for k = 1:N_FREQ_POINTS
        z_k = exp(1i * frequency_vector(k));

        Tzw_l1     = evalfr_FIR(lft_l1_jth_channel, z_k);
        Tzw_spreg  = evalfr_FIR(lft_spreg_jth_channel, z_k);
        Tzw_oracle = evalfr_FIR(lft_oracle_jth_channel, z_k);

        cost_l1(k)     = norm(Tzw_l1, Inf);
        cost_spreg(k)  = norm(Tzw_spreg, Inf);
        cost_oracle(k) = norm(Tzw_oracle, Inf);

        % assert(cost_l1(k) == norm(Tzw_l1, Inf), "L1 cost is not equal to the infinity norm " + num2str(k));
    end
end
% Before plotting, we print the sum of the costs for each controller
disp("Plotting channel " + num2str(input_channel) + ".")
% fprintf("Channel %d:  --> L1 cost = %.4f, SpReg cost = %.4f, Oracle cost = %.4f", input_channel, sum(cost_l1), sum(cost_spreg), sum(cost_oracle));
% fprintf(" (L1:  %.2f%% wrt Spreg)\n", 100 * (sum(cost_l1) - sum(cost_spreg)) / sum(cost_spreg));


color_l1 = plot_params.color_l1;
color_spreg = plot_params.color_spreg;
color_oracle = plot_params.color_oracle;

linestyle_l1 = plot_params.linestyle_l1;
linestyle_spreg = plot_params.linestyle_spreg;
linestyle_oracle = plot_params.linestyle_oracle;

% line_width = plot_params.line_width;
with_logy = plot_params.with_logy;

final_figure = figure;
hold on;

if with_logy
    plot_l1 = plot(frequency_vector, log10(cost_l1), 'Color', color_l1, 'LineStyle', linestyle_l1, 'LineWidth', plot_params.line_width_l1, 'DisplayName', "$\mathbf{K}^{\mathcal{L}_{1}}$");
    plot_spreg = plot(frequency_vector, log10(cost_spreg), 'Color', color_spreg, 'LineStyle', linestyle_spreg, 'LineWidth', plot_params.line_width_spreg, 'DisplayName', "$\mathbf{K}^{{SR}}_{\mathcal{L}_{1}}$");
    if plot_params.plot_the_oracle        
        plot_oracle = plot(frequency_vector, log10(cost_oracle), 'Color', color_oracle, 'LineStyle', linestyle_oracle, 'LineWidth', plot_params.line_width_oracle, 'DisplayName', "$\mathbf{\hat{K}}$");
    end
else
    plot_l1 = plot(frequency_vector, cost_l1, 'Color', color_l1, 'LineStyle', linestyle_l1, 'LineWidth', plot_params.line_width_l1, 'DisplayName', "$\mathbf{K}^{\mathcal{L}_{1}}$");
    plot_spreg = plot(frequency_vector, cost_spreg, 'Color', color_spreg, 'LineStyle', linestyle_spreg, 'LineWidth', plot_params.line_width_spreg, 'DisplayName', "$\mathbf{K}^{{SR}}_{\mathcal{L}_{1}}$");
    if plot_params.plot_the_oracle        
        plot_oracle = plot(frequency_vector, cost_oracle, 'Color', color_oracle, 'LineStyle', linestyle_oracle, 'LineWidth', plot_params.line_width_oracle, 'DisplayName', "$\mathbf{\hat{K}}$");
    end
end
if plot_params.plot_the_oracle
    legend([plot_spreg, plot_l1, plot_oracle], 'Interpreter','latex', 'FontSize',16)
else
    legend([plot_spreg, plot_l1], 'Interpreter','latex', 'FontSize',16);
end
if with_logy
    ylabel('$\log_{10}(||\mathbf{F}^{[:,' + string(input_channel) + ']}(e^{j\Omega})||)_{\infty}^2$', 'Interpreter','latex', 'FontSize',16);
else
    ylabel('$||\mathbf{F}^{[:,' + string(input_channel) + ']}(e^{j\Omega})||_{\infty}^2$', 'Interpreter','latex', 'FontSize',16);
end
grid minor;
if half_frequencies
    xlim([0, pi]);
    xticks([0, pi/4, pi/2, 3*pi/4, pi]);
    xticklabels({'0', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3 \pi}{4}$', '$\pi$'});
else
    xlim([-pi, pi]);
    xticks([-pi, -pi/2, 0, pi/2, pi]);
    xticklabels({'$-\pi$', '$-\frac{\pi}{2}$', '0', '$\frac{\pi}{2}$', '$\pi$'});
end
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
xaxisproperties.FontSize = 12;
xlabel('Frequency $\omega$', 'Interpreter','latex', 'FontSize', 16);

if ~isempty(plot_params.Position)
    % Set the figure position if specified
    final_figure.Position = plot_params.Position;
end
hold off;
end



function Phis_mat_3D = compute_FIR_Matrix(sys, Matrices_struct, N_tf)
% Computes the FIR transfer matrix Phis_L1_mat from sys and Phis_l1.mat

n = sys.n;
p = sys.p;

Phis_mat_3D = zeros(sys.n_z, sys.n_w, N_tf + 1);

Rmat = Matrices_struct.R;
Mmat = Matrices_struct.M;
Nmat = Matrices_struct.N;
Lmat = Matrices_struct.L;

Tzw_0 = sys.D12 * Lmat(:,1:p) * sys.D21 + sys.D11;
Phis_mat_3D(:,:,1) = Tzw_0;

for t = 1:N_tf
    Tzw_t = sys.C1  * Rmat(:, t*n+1:(t+1)*n) * sys.B1 + ...
            sys.D12 * Mmat(:, t*n+1:(t+1)*n) * sys.B1 + ...
            sys.C1  * Nmat(:, t*p+1:(t+1)*p) * sys.D21 + ...
            sys.D12 * Lmat(:, t*p+1:(t+1)*p) * sys.D21;
    Phis_mat_3D(:,:,t+1) = Tzw_t;
end
end



function [evaluated_matrix] = evalfr_FIR(Phi_FIR , z_k)

[~,~, N_tf_plus_one ] = size(Phi_FIR);
evaluated_matrix = zeros(size(Phi_FIR,1), size(Phi_FIR,2));
for i = 1:N_tf_plus_one
    evaluated_matrix = evaluated_matrix + Phi_FIR(:,:,i) * (z_k)^(i-1);
end
end