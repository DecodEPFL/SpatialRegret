%% plots_for_SDP.m
% Generates frequency response plots for SDP-based controller comparison
%
% This script creates frequency domain plots comparing different controllers
%
% The script generates plots for each input channel and a combined plot.
% Figures are saved in the specified folder as both PDF and FIG files.
%
% Required variables (must be defined before running this script):
%   - lft_oracle, lft_spreg, lft_h2, lft_hinf: Closed-loop transfer functions
%   - n_agents: Number of agents in the network
%   - savingfolder_figures: Directory path for saving figures
%   - chosen_spreg: Index of the spatial regret controller to plot
%
% See also: plots_for_L1, draw_frequency_plots

%% Configuration
% Construct the plot parameters structure
plot_params = struct();
plot_params.N_FREQ_POINTS = 1000;       % Number of frequency points
plot_params.input_channel = 0;          % Input channel to plot (0 for all)
plot_params.half_frequencies = true;    % Use half frequencies (0 to pi)

% Color scheme for different controllers
plot_params.color_h2     = [0.2157, 0.4941, 0.7225];  % Blue
plot_params.color_hinf   = [0.9290, 0.6940, 0.1250];  % Yellow
plot_params.color_spreg  = [0.8941, 0.1019, 0.1098];  % Red
plot_params.color_oracle = 'k';                       % Black

% Line styles
plot_params.linestyle_h2     = '-';
plot_params.linestyle_hinf   = '-';
plot_params.linestyle_spreg  = '-';
plot_params.linestyle_oracle = '--';

% Line widths
line_width = 2.5;
plot_params.line_width_h2     = line_width;
plot_params.line_width_hinf   = line_width;
plot_params.line_width_spreg  = line_width;
plot_params.line_width_oracle = 0.6 * line_width;

% Plot settings
plot_params.with_logy = false;                 % Use logarithmic scale for y-axis
plot_params.Position  = [100, 100, 556, 272];  % Figure position and size

% Check if save_plots_as_figs is defined, default to true
if ~exist('save_plots_as_figs', 'var')
    save_plots_as_figs = true;
end

% Create output directory if it doesn't exist
if ~exist(savingfolder_figures, 'dir')
    mkdir(savingfolder_figures)
end

%% Generate plots for each individual channel
fprintf('Generating frequency response plots...\n');

for input_channel = 1:n_agents
    fprintf('  Output Channel %d/%d\n', input_channel, n_agents);
    
    % Set the input channel for the current iteration
    plot_params.input_channel = input_channel;
    
    % Generate frequency plot
    spreg_spectr_norm_Fig_i = draw_frequency_plots(lft_spreg, lft_oracle, ...
                                                    lft_h2, lft_hinf, plot_params);
    
    % Save as PDF (vector format)
    filename_base = sprintf('Spectral_freq_plot_on_channel_%d', input_channel);
    exportgraphics(spreg_spectr_norm_Fig_i, ...
                   savingfolder_figures + filename_base + ".pdf", ...
                   'ContentType', 'vector');
    
    % Save as FIG (MATLAB format) if requested
    if save_plots_as_figs
        saveas(spreg_spectr_norm_Fig_i, ...
               savingfolder_figures + filename_base + ".fig");
    end
end

%% Generate combined plot for all channels
fprintf('  All channels (combined)\n');

% Set to 0 to plot all channels
plot_params.input_channel = 0;

% Generate frequency plot
spreg_spectr_norm_Fig_all = draw_frequency_plots(lft_spreg, lft_oracle, ...
                                                  lft_h2, lft_hinf, plot_params);

% Save as PDF
exportgraphics(spreg_spectr_norm_Fig_all, ...
               savingfolder_figures + "Spectral_freq_plot_all_channels.pdf", ...
               'ContentType', 'vector');

% Save as FIG if requested
if save_plots_as_figs
    saveas(spreg_spectr_norm_Fig_all, ...
           savingfolder_figures + "Spectral_freq_plot_all_channels.fig");
end

fprintf('Frequency response plots saved to: %s\n', savingfolder_figures);


%% Helper Function: draw_frequency_plots
function [spreg_spectr_norm_Fig] = draw_frequency_plots(lft_spreg, lft_oracle, ...
                                                         lft_h2, lft_hinf, plot_params)
% DRAW_FREQUENCY_PLOTS - Generate frequency response comparison plots
%
% Inputs:
%   lft_spreg  - Spatial regret closed-loop transfer function
%   lft_oracle - Oracle closed-loop transfer function
%   lft_h2     - H2-optimal closed-loop transfer function
%   lft_hinf   - H-infinity closed-loop transfer function
%   plot_params - Structure with plotting parameters
%
% Output:
%   spreg_spectr_norm_Fig - Figure handle

% Extract parameters
half_frequencies = plot_params.half_frequencies;    
N_FREQ_POINTS    = plot_params.N_FREQ_POINTS;
input_channel    = plot_params.input_channel;

% Define frequency vector
if half_frequencies
    frequency_vector = linspace(0, pi, N_FREQ_POINTS);
else
    frequency_vector = linspace(-pi, pi, N_FREQ_POINTS);
end

% Initialize cost arrays for the different controllers
cost_h2     = zeros(1, N_FREQ_POINTS);
cost_hinf   = zeros(1, N_FREQ_POINTS);
cost_spreg  = zeros(1, N_FREQ_POINTS);
cost_oracle = zeros(1, N_FREQ_POINTS);

% Extract specific input channel if requested (0 means all channels)
if input_channel > 0 
    lft_h2     = lft_h2(:, input_channel);
    lft_hinf   = lft_hinf(:, input_channel);
    lft_spreg  = lft_spreg(:, input_channel);
    lft_oracle = lft_oracle(:, input_channel);
end

% Compute frequency response at each frequency point
for k = 1:N_FREQ_POINTS
    z_k = exp(1i * frequency_vector(k));
    
    % Evaluate transfer functions at frequency z_k
    Tzw_h2     = evalfr(lft_h2, z_k);
    Tzw_hinf   = evalfr(lft_hinf, z_k);
    Tzw_spreg  = evalfr(lft_spreg, z_k);
    Tzw_oracle = evalfr(lft_oracle, z_k);
    
    % Compute singular values
    svd_h2     = svd(Tzw_h2);
    svd_hinf   = svd(Tzw_hinf);
    svd_spreg  = svd(Tzw_spreg);
    svd_oracle = svd(Tzw_oracle);
    
    % Store squared maximum singular value
    cost_h2(k)     = svd_h2(1)^2;
    cost_hinf(k)   = svd_hinf(1)^2;
    cost_spreg(k)  = svd_spreg(1)^2;
    cost_oracle(k) = svd_oracle(1)^2;
end

% Extract plotting parameters
color_h2     = plot_params.color_h2;
color_hinf   = plot_params.color_hinf;
color_spreg  = plot_params.color_spreg;
color_oracle = plot_params.color_oracle;

linestyle_h2     = plot_params.linestyle_h2;
linestyle_hinf   = plot_params.linestyle_hinf;
linestyle_spreg  = plot_params.linestyle_spreg;
linestyle_oracle = plot_params.linestyle_oracle;

line_width_h2     = plot_params.line_width_h2;
line_width_hinf   = plot_params.line_width_hinf;
line_width_spreg  = plot_params.line_width_spreg;
line_width_oracle = plot_params.line_width_oracle;

with_logy = plot_params.with_logy;

% Create figure and plot frequency responses
spreg_spectr_norm_Fig = figure;
hold on;

if with_logy
    % Logarithmic scale
    plot_h2 = plot(frequency_vector, log10(cost_h2), ...
                   'Color', color_h2, 'LineStyle', linestyle_h2, ...
                   'LineWidth', line_width_h2, ...
                   'DisplayName', "$\mathbf{K}^{\mathcal{H}_{2}}$");
    plot_hinf = plot(frequency_vector, log10(cost_hinf), ...
                     'Color', color_hinf, 'LineStyle', linestyle_hinf, ...
                     'LineWidth', line_width_hinf, ...
                     'DisplayName', "$\mathbf{K}^{\mathcal{H}_{\infty}}$");
    plot_oracle = plot(frequency_vector, log10(cost_oracle), ...
                       'Color', color_oracle, 'LineStyle', linestyle_oracle, ...
                       'LineWidth', line_width_oracle, ...
                       'DisplayName', "$\mathbf{\hat{K}}$");
    plot_spreg = plot(frequency_vector, log10(cost_spreg), ...
                      'Color', color_spreg, 'LineStyle', linestyle_spreg, ...
                      'LineWidth', line_width_spreg, ...
                      'DisplayName', "$\mathbf{K}^{{SR}}$");
else
    % Linear scale
    plot_h2 = plot(frequency_vector, cost_h2, ...
                   'Color', color_h2, 'LineStyle', linestyle_h2, ...
                   'LineWidth', line_width_h2, ...
                   'DisplayName', "$\mathbf{K}^{\mathcal{H}_{2}}$");
    plot_hinf = plot(frequency_vector, cost_hinf, ...
                     'Color', color_hinf, 'LineStyle', linestyle_hinf, ...
                     'LineWidth', line_width_hinf, ...
                     'DisplayName', "$\mathbf{K}^{\mathcal{H}_{\infty}}$");
    plot_oracle = plot(frequency_vector, cost_oracle, ...
                       'Color', color_oracle, 'LineStyle', linestyle_oracle, ...
                       'LineWidth', line_width_oracle, ...
                       'DisplayName', "$\mathbf{\hat{K}}$");
    plot_spreg = plot(frequency_vector, cost_spreg, ...
                      'Color', color_spreg, 'LineStyle', linestyle_spreg, ...
                      'LineWidth', line_width_spreg, ...
                      'DisplayName', "$\mathbf{K}^{{SR}}$");
end

% Add legend with custom order (Spatial Regret, H2, Hinf, Oracle)
legend([plot_spreg, plot_h2, plot_hinf, plot_oracle], ...
       'Interpreter', 'latex', 'FontSize', 16);

% Set y-axis label based on channel and scale
if input_channel > 0
    % Specific channel
    if with_logy
        ylabel('$\log_{10}(||\mathbf{F}^{[:,' + string(input_channel) + ']}(e^{j\omega})||_2^2)$', 'Interpreter','latex', 'FontSize',18);
    else
        ylabel('$||\mathbf{F}^{[:,' + string(input_channel) + ']}(e^{j\omega})||_2^2$', 'Interpreter','latex', 'FontSize',18);
    end
else
    % All channels
    if with_logy
        ylabel('$\log_{10}(||\mathbf{F}_{\ell}(e^{j\omega})||_2^2)$', 'Interpreter','latex', 'FontSize',16);
    else
        ylabel('$||\mathbf{F}_{\ell}(e^{j\omega})||_2^2$', 'Interpreter','latex', 'FontSize',16);
    end
end

% Grid settings
grid minor;

% Set x-axis limits and ticks
if half_frequencies
    xlim([0, pi]);
    xticks([0, pi/4, pi/2, 3*pi/4, pi]);
    xticklabels({'0', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', ...
                 '$\frac{3 \pi}{4}$', '$\pi$'});
else
    xlim([-pi, pi]);
    xticks([-pi, -pi/2, 0, pi/2, pi]);
    xticklabels({'$-\pi$', '$-\frac{\pi}{2}$', '0', ...
                 '$\frac{\pi}{2}$', '$\pi$'});
end

% Configure x-axis properties
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
xaxisproperties.FontSize = 12;

% Set x-axis label
xlabel('Frequency $\omega$', 'Interpreter', 'latex', 'FontSize', 16);

% Set figure position and size
spreg_spectr_norm_Fig.Position = plot_params.Position;

hold off;
end