%% Sim_one_specific_disturbance_SDP.m
% Simulates time-domain response to a specific disturbance for SDP-based controllers
%
% This script generates a manual disturbance signal at specified frequencies
% and simulates the closed-loop response for different controllers:
%   - H2-optimal controller
%   - H-infinity controller
%   - Oracle controller
%   - Spatial Regret controller
%
% The script computes and compares the performance costs and generates
% time-domain plots showing the disturbance and instantaneous costs.
%
% Required variables (must be defined before running this script):
%   - lft_h2, lft_hinf, lft_oracle, lft_spreg: Closed-loop transfer functions
%   - sys: System structure with plant model
%   - chosen_channel: Channel index for disturbance injection
%   - savingfolder_figures: Directory path for saving figures
%   - Ts: Sampling time
%
% See also: Sim_one_specific_disturbance_L1, plots_for_SDP

%% Plot Configuration
% Transparency settings
alpha_base  = 1.0;
alpha_spreg = 1.0;

% Plot parameters structure
plot_params = struct;
plot_params.line_width_base  = 2;
plot_params.line_width_spreg = 2;
plot_params.with_logy = false;
plot_params.Position  = [642, 337, 372, 310];

% Line styles
plot_params.linestyle_h2     = '-';
plot_params.linestyle_hinf   = '-';
plot_params.linestyle_spreg  = '-';
plot_params.linestyle_oracle = '--';

% Color scheme
plot_params.color_h2     = [0.2157, 0.4941, 0.7225];  % Blue
plot_params.color_hinf   = [0.9290, 0.6940, 0.1250];  % Yellow
plot_params.color_spreg  = [0.8941, 0.1019, 0.1098];  % Red
plot_params.color_oracle = [0.4660, 0.6740, 0.1880];  % Green

% Check if save_plots_as_figs is defined, default to true
if ~exist('save_plots_as_figs', 'var')
    save_plots_as_figs = true;
end

% Create output directory if it doesn't exist
if ~exist(savingfolder_figures, 'dir')
    mkdir(savingfolder_figures)
end

%% Generate Disturbance Signal
% Define disturbance frequencies (normalized by sampling time)
desired_normalized_omegas = [0.1, pi/5];  % Example: multiple frequency components

% Simulation time parameters
T_END = 25;                          % End time (seconds)
t = 0:sys.plant.Ts:T_END;            % Time vector
n_instants = length(t);              % Number of time steps

% Initialize disturbance signal
w_manual = zeros(n_instants, sys.n_w);

% Generate multi-frequency disturbance
for freq = desired_normalized_omegas
    omega_desired = freq / Ts;
    
    % Scale DC amplitude to be 1/2 of other frequencies
    if freq == 0
        A = 1/2;
    else
        A = 1;
    end
    
    % Add cosine component at specified frequency
    w_manual(:, chosen_channel) = w_manual(:, chosen_channel) + ...
                                  A * cos(omega_desired * t)';
end

% Normalize the disturbance signal
w_manual = w_manual / norm(w_manual);

%% Simulate Closed-Loop Responses
fprintf('Simulating closed-loop responses...\n');

z_h2     = lsim(lft_h2, w_manual, t);
z_hinf   = lsim(lft_hinf, w_manual, t);
z_oracle = lsim(lft_oracle, w_manual, t);
z_spreg  = lsim(lft_spreg, w_manual, t);

fprintf('Simulation complete.\n');

%% Compute Performance Costs
% Total cost (sum of squared outputs)
cost_h2     = sum(sum(z_h2.^2));
cost_hinf   = sum(sum(z_hinf.^2));
cost_oracle = sum(sum(z_oracle.^2));
cost_spreg  = sum(sum(z_spreg.^2));

% Instantaneous cost (sum of squared outputs at each time)
cost_h2_time     = sum(z_h2.^2, 2);
cost_hinf_time   = sum(z_hinf.^2, 2);
cost_oracle_time = sum(z_oracle.^2, 2);
cost_spreg_time  = sum(z_spreg.^2, 2);

% Display results
fprintf('\n================================================\n');
fprintf('PERFORMANCE COMPARISON - Total Costs\n');
fprintf('================================================\n');
fprintf('Cost H2:     %.4f\n', cost_h2);
fprintf('Cost Hinf:   %.4f\n', cost_hinf);
fprintf('Cost Oracle: %.4f\n', cost_oracle);
fprintf('Cost SpReg:  %.4f\n', cost_spreg);
fprintf('================================================\n\n');

%% Plot 1: Disturbance and Instantaneous Costs
figure;

% Subplot 1: Disturbance signal
subplot(2, 1, 1);
plot(t, w_manual(:, chosen_channel), 'LineWidth', 2, ...
     'DisplayName', "$w_t^{[" + num2str(chosen_channel) + "]}$");
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Disturbance', 'Interpreter', 'latex');
title("Manual Disturbance on Agent " + num2str(chosen_channel));
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 16);
grid on;

% Subplot 2: Instantaneous costs
subplot(2, 1, 2);
p_hinf = plot(t, cost_hinf_time, ...
              'LineWidth', plot_params.line_width_base, ...
              'DisplayName', "$\mathbf{K}^{{H}_{\infty}}$", ...
              'Color', plot_params.color_hinf);
hold on;
p_h2 = plot(t, cost_h2_time, ...
            'LineWidth', plot_params.line_width_base, ...
            'DisplayName', "$\mathbf{K}^{{H}_{2}}$", ...
            'Color', plot_params.color_h2);
p_spreg = plot(t, cost_spreg_time, ...
               'LineWidth', plot_params.line_width_spreg, ...
               'DisplayName', "$\mathbf{K}^{{SR}}$", ...
               'Color', plot_params.color_spreg);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Cost', 'Interpreter', 'latex');
legend([p_spreg, p_h2, p_hinf], 'Location', 'best', ...
       'Interpreter', 'latex', 'FontSize', 16);
grid on;
hold off;

% Overall title
sgtitle("DISTURBANCE WITH FREQs [ " + num2str(desired_normalized_omegas) + " ]");

% Save figure
exportgraphics(gcf, savingfolder_figures + ...
               "manual_disturbance_on_mass_" + num2str(chosen_channel) + ".pdf", ...
               'ContentType', 'vector');
if save_plots_as_figs
    saveas(gcf, savingfolder_figures + ...
           "manual_disturbance_on_mass_" + num2str(chosen_channel) + ".fig");
end

%% Compute and Display Mean Costs
mean_cost_h2     = mean(cost_h2_time);
mean_cost_hinf   = mean(cost_hinf_time);
mean_cost_spreg  = mean(cost_spreg_time);
mean_cost_oracle = mean(cost_oracle_time);

% fprintf('Mean Costs:\n');
% fprintf('  H2:     %.4f\n', mean_cost_h2);
% fprintf('  Hinf:   %.4f\n', mean_cost_hinf);
% fprintf('  SpReg:  %.4f\n', mean_cost_spreg);
% fprintf('  Oracle: %.4f\n', mean_cost_oracle);
% fprintf('\n');

% Compute improvement percentages
improv_over_h2 = (mean_cost_h2 - mean_cost_spreg) / mean_cost_h2 * 100;
improv_over_hinf = (mean_cost_hinf - mean_cost_spreg) / mean_cost_hinf * 100;

fprintf('Spatial Regret Improvements:\n');
fprintf('  Over H2:   %.2f%%\n', improv_over_h2);
fprintf('  Over Hinf: %.2f%%\n', improv_over_hinf);
fprintf('\n');

%% Plot 2: Cost Comparison Only (Publication Quality)
cost_manual_disturbance_fig = figure;

p_hinf = plot(t, cost_hinf_time, ...
              'LineWidth', plot_params.line_width_base, ...
              'DisplayName', "$\mathbf{K}^{{H}_{\infty}}$", ...
              'Color', [plot_params.color_hinf, alpha_base]);
hold on;
p_h2 = plot(t, cost_h2_time, ...
            'LineWidth', plot_params.line_width_base, ...
            'DisplayName', "$\mathbf{K}^{{H}_{2}}$", ...
            'Color', [plot_params.color_h2, alpha_base]);
p_spreg = plot(t, cost_spreg_time, ...
               'LineWidth', plot_params.line_width_spreg, ...
               'DisplayName', "$\mathbf{K}^{{SR}}$", ...
               'Color', [plot_params.color_spreg, alpha_spreg]);

xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$||z_t||^2$', 'Interpreter', 'latex', 'FontSize', 18);
legend([p_spreg, p_h2, p_hinf], 'Location', 'northeast', ...
       'Interpreter', 'latex', 'FontSize', 16);

cost_manual_disturbance_fig.Position = [714, 494, 590, 186];
grid on;
hold off;

% Save figure
exportgraphics(gcf, savingfolder_figures + ...
               "only_cost_manual_disturbance_on_mass_" + num2str(chosen_channel) + ".pdf", ...
               'ContentType', 'vector');
if save_plots_as_figs
    saveas(gcf, savingfolder_figures + ...
           "only_cost_manual_disturbance_on_mass_" + num2str(chosen_channel) + ".fig");
end

fprintf('Plots saved to: %s\n', savingfolder_figures);

%% Optional: Plot State and Input Trajectories (Currently Commented Out)
% %% Saving also the plots of the states
% saving_folder_plots_states = savingfolder_figures + "/States_trajectories/";
% if ~exist(saving_folder_plots_states, 'dir')
%     mkdir(saving_folder_plots_states)
% end
% 
% 
% 
% % First we extract the states from the z vectors
% 
% states_indeces  = 1:2*n_agents;
% input_indeces = 2*n_agents + 1:3*n_agents;
% 
% states_h2 = z_h2(:, states_indeces);
% states_hinf = z_hinf(:, states_indeces);
% states_oracle = z_oracle(:, states_indeces);
% states_spreg = z_spreg(:, states_indeces);
% 
% inputs_h2 = z_h2(:, input_indeces);
% inputs_hinf = z_hinf(:, input_indeces);
% inputs_oracle = z_oracle(:, input_indeces);
% inputs_spreg = z_spreg(:, input_indeces);
% 
% % Plot the states trajectories
% figure;
% subplot(2,1,1);
% plot(t, states_h2, 'LineWidth', plot_params.line_width_base, 'DisplayName', "$\mathbf{K}^{{H}_{2}}$", 'Color', plot_params.color_h2);
% hold on;
% plot(t, states_hinf, 'LineWidth', plot_params.line_width_base, 'DisplayName', "$\mathbf{K}^{{H}_{\infty}}$", 'Color', plot_params.color_hinf);
% plot(t, states_spreg, 'LineWidth', plot_params.line_width_spreg, 'DisplayName', "$\mathbf{K}^{{SR}}$", 'Color', plot_params.color_spreg);
% xlabel('Time [s]');
% ylabel('States');
% title("States trajectories for disturbance on mass " + num2str(chosen_channel));
% % legend('Location','best', 'Interpreter','latex', 'FontSize',16);
% 
% subplot(2,1,2);
% plot(t, inputs_h2, 'LineWidth', plot_params.line_width_base, 'DisplayName', "$\mathbf{K}^{{H}_{2}}$", 'Color', plot_params.color_h2);
% hold on;
% plot(t, inputs_hinf, 'LineWidth', plot_params.line_width_base, 'DisplayName', "$\mathbf{K}^{{H}_{\infty}}$", 'Color', plot_params.color_hinf);
% plot(t, inputs_spreg, 'LineWidth', plot_params.line_width_spreg, 'DisplayName', "$\mathbf{K}^{{SR}}$", 'Color', plot_params.color_spreg);
% xlabel('Time [s]');
% ylabel('Inputs');
% title("Inputs trajectories for disturbance on mass " + num2str(chosen_channel));
% % legend('Location','best', 'Interpreter','latex', 'FontSize',16);
% grid on;
% hold off;
% sgtitle("States and Inputs trajectories for disturbance with Freqs [ " + num2str(desired_normalized_omegas) + " ]");
% exportgraphics(gcf, saving_folder_plots_states + "states_inputs_trajectories_on_mass_" + num2str(chosen_channel) + ".pdf", 'ContentType', 'vector');
% saveas(gcf, saving_folder_plots_states + "states_inputs_trajectories_on_mass_" + num2str(chosen_channel) + ".fig");
% 
% % Plot the states trajectories for each agent
% for agent = 1:n_agents
%     figure;
%     subplot(2,1,1);
%     plot(t, states_h2(:, (agent-1)*2+1), 'LineWidth', plot_params.line_width_base, 'DisplayName', "$\mathbf{K}^{{H}_{2}}$ - Agent " + num2str(agent), 'Color', plot_params.color_h2);
%     hold on;
%     plot(t, states_hinf(:, (agent-1)*2+1), 'LineWidth', plot_params.line_width_base, 'DisplayName', "$\mathbf{K}^{{H}_{\infty}}$ - Agent " + num2str(agent), 'Color', plot_params.color_hinf);
%     plot(t, states_spreg(:, (agent-1)*2+1), 'LineWidth', plot_params.line_width_spreg, 'DisplayName', "$\mathbf{K}^{{SR}}$ - Agent " + num2str(agent), 'Color', plot_params.color_spreg);
%     xlabel('Time [s]');
%     ylabel('Position');
%     title("Position trajectory for disturbance on mass " + num2str(chosen_channel) + " - Agent " + num2str(agent));
%     legend('Location','best', 'Interpreter','latex', 'FontSize',16);
% 
%     subplot(2,1,2);
%     plot(t, states_h2(:, (agent-1)*2+2), 'LineWidth', plot_params.line_width_base, 'DisplayName', "$\mathbf{K}^{{H}_{2}}$ - Agent " + num2str(agent), 'Color', plot_params.color_h2);
%     hold on;
%     plot(t, states_hinf(:, (agent-1)*2+2), 'LineWidth', plot_params.line_width_base, 'DisplayName', "$\mathbf{K}^{{H}_{\infty}}$ - Agent " + num2str(agent), 'Color', plot_params.color_hinf);
%     plot(t, states_spreg(:, (agent-1)*2+2), 'LineWidth', plot_params.line_width_spreg, 'DisplayName', "$\mathbf{K}^{{SR}}$ - Agent " + num2str(agent), 'Color', plot_params.color_spreg);
%     xlabel('Time [s]');
%     ylabel('Velocity');
%     title("Velocity trajectory for disturbance on mass " + num2str(chosen_channel)+ " - Agent " + num2str(agent));
%     legend('Location','best', 'Interpreter','latex', 'FontSize',16);
%     grid on;
%     hold off;
%     sgtitle("States trajectories for disturbance with Freqs [ " + num2str(desired_normalized_omegas) + " ] - Agent " + num2str(agent));
%     exportgraphics(gcf, saving_folder_plots_states + "states_trajectories_on_mass_" + num2str(chosen_channel) + "_agent_" + num2str(agent) + ".pdf", 'ContentType', 'vector');
%     saveas(gcf, saving_folder_plots_states + "states_trajectories_on_mass_" + num2str(chosen_channel) + "_agent_" + num2str(agent) + ".fig");
% end
