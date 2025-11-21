%Parameters for the plot

% Define the different parameter for the plots
alpha_L1 = 1.0;
alpha_Spreg_L1 = 1.0;
plot_params = struct;

plot_params.N_FREQ_POINTS = 200;
plot_params.line_width_l1 = 2;
plot_params.line_width_spreg = 2;
% plot_params.line_width_oracle = 2;

plot_params.with_logy = false;
plot_params.Position = [642,337,372,310];

plot_params.linestyle_l1 = '-';
plot_params.linestyle_spreg = '-';
% plot_params.linestyle_oracle = '-.';

plot_params.color_l1     = [0.2157, 0.4941, 0.7225];   % blue
plot_params.color_spreg  = [0.8941, 0.1019, 0.1098];   % red


Phis_L1_mat    = compute_FIR_Matrix(sys, Phis_l1.mat, N_tf);
Phis_Spreg_mat = compute_FIR_Matrix(sys, Phis_spreg.mat, N_tf);

if ~exist(savingfolder_figures, 'dir')
    mkdir(savingfolder_figures)
end

desired_normalized_omegas = linspace(3/4*pi, pi, 100); disp("100 frequencies from 3/4*pi to pi");

T_END = 40;
t = 0:sys.plant.Ts:T_END;
n_instants = length(t);
N_EXPERIMENTS = 1e4;

improvement_vector = zeros(1,N_EXPERIMENTS);

parfor i = 1:N_EXPERIMENTS 
w_manual = zeros(n_instants, sys.n_w);

for freq = desired_normalized_omegas
    omega_desired = freq / Ts;
    phase = 2*pi*rand(); % random phase in [0, 2pi]
    if freq == 0
        A = 1/2;
    else
        A = 1;
    end
    w_manual(:, chosen_channel) = w_manual(:, chosen_channel) + A * cos(omega_desired * t' + phase);
end

w_manual = w_manual / max(abs(w_manual), [], 'all');

% z_l1 = lsim(lft_l1, w_manual, t);
% z_spreg_l1 = lsim(lft_spreg, w_manual, t);
z_l1       = simulateFIR(Phis_L1_mat, w_manual);
z_spreg_l1 = simulateFIR(Phis_Spreg_mat, w_manual);

cost_l1       = max(abs(z_l1),[],2)/max(w_manual, [], "all");
cost_spreg_l1 = max(abs(z_spreg_l1),[],2)/max(w_manual, [], "all");


l1_cost_l1 = max(cost_l1);
l1_cost_l1_spreg = max(cost_spreg_l1);

ith_improvement =(l1_cost_l1 - l1_cost_l1_spreg)/l1_cost_l1*100;
improvement_vector(i) = ith_improvement;
end


disp("Average improvement over L1: " + num2str(mean(improvement_vector)) + "%");

%% Plots
w_manual = zeros(n_instants, sys.n_w);

for freq = desired_normalized_omegas
    omega_desired = freq / Ts;
    phase = 2*pi*rand(); % random phase in [0, 2pi]
    if freq == 0
        A = 1/2;
    else
        A = 1;
    end
    w_manual(:, chosen_channel) = w_manual(:, chosen_channel) + A * cos(omega_desired * t' + phase);
end

z_l1       = simulateFIR(Phis_L1_mat, w_manual);
z_spreg_l1 = simulateFIR(Phis_Spreg_mat, w_manual);

cost_l1       = max(abs(z_l1),[],2)/max(w_manual, [], "all");
cost_spreg_l1 = max(abs(z_spreg_l1),[],2)/max(w_manual, [], "all");


l1_cost_l1 = max(cost_l1);
l1_cost_l1_spreg = max(cost_spreg_l1);


figure;
subplot(2,1,1);
plot(t, w_manual(:,chosen_channel), 'LineWidth', 2, 'DisplayName', "$w_t^{[" + num2str(chosen_channel) + "]}$");
xlabel('Time [s]');
ylabel('Disturbance');
title("Manual Disturbance on the mass no. " + num2str(chosen_channel));
legend('Location','best', 'Interpreter','latex', 'FontSize',16);
grid on;

subplot(2,1,2);
plot_l1 = plot(t, cost_l1, 'LineWidth', plot_params.line_width_l1, 'DisplayName', "$\mathbf{K}^{{L}_{1}}$", 'Color', plot_params.color_l1);
hold on;
plot_spreg = plot(t, cost_spreg_l1, 'LineWidth', plot_params.line_width_spreg, 'DisplayName', "$\mathbf{K}^{{SR}}$", 'Color', plot_params.color_spreg);
xlabel('Time [s]');
ylabel('Cost');
title("Instantaneous costs (Improv over L1: " + num2str((l1_cost_l1 - l1_cost_l1_spreg)/l1_cost_l1*100) + "%)");
legend([plot_spreg, plot_l1], 'Location','best','Interpreter','latex', 'FontSize',16);
grid on;
hold off;

% sgtitle("DISTURBANCE WITH FREQs [ " + num2str(desired_normalized_omegas) + " ]");
exportgraphics(gcf, savingfolder_figures + "manual_disturbance_on_mass_" + num2str(chosen_channel) + ".pdf", 'ContentType', 'vector');



cost_manual_disturbance_fig = figure;
plot_l1 = plot(t, cost_l1, 'LineWidth', plot_params.line_width_l1, 'DisplayName', "$\mathbf{K}^{{L}_{1}}$", 'Color', [plot_params.color_l1, alpha_L1]);
hold on;
plot_spreg = plot(t, cost_spreg_l1, 'LineWidth', plot_params.line_width_spreg, 'DisplayName', "$\mathbf{K}^{{SR}}$", 'Color', [plot_params.color_spreg, alpha_Spreg_L1] );
xlabel('Time $t$ [s]', 'Interpreter','latex', 'FontSize',18);
ylabel('$||z_t||_{\infty}$', 'Interpreter','latex', 'FontSize',18);
legend([plot_spreg, plot_l1], 'Location','best','Interpreter','latex', 'FontSize',16);
cost_manual_disturbance_fig.Position = [714   494   590   186];
grid on;
hold off;
exportgraphics(gcf, savingfolder_figures + "only_cost_manual_disturbance_on_mass_" + num2str(chosen_channel) + ".pdf", 'ContentType', 'vector');
saveas(gcf, savingfolder_figures + "only_cost_manual_disturbance_on_mass_" + num2str(chosen_channel) + ".fig");







%% AUXILIARY LOCAL FUNCTIONS


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





