function [K_opt, Q_opt, objective] = calculus_distributed_H2_youla_sampled_WITH_TF(sys, delays_matrix,options)

P11 = sys.P11;
P12 = sys.P12;
P21 = sys.P21;

%% Optimization Variables

% Q_sdp = sdpvar(sys.m,sys.p,options.N_tf+1);
Q_sdp = cell(options.N_tf+1,1);

%% Information Sparsities
N_agents = size(delays_matrix,1);
m_i = sys.m/N_agents;
p_i = sys.p/N_agents;
Q_mask = ones(sys.m,sys.p,options.N_tf+1);

for i=1:N_agents        %rows --> p_i
    for j=1:N_agents    %columns --> m_i
        dist_ij = delays_matrix(i,j);
        %HERE, if the sub-systems are not homogeneous, you need to compute
        % m_i = Partition_inputs_matrix(i,j)
        % p_i = Partition_outputs_matrix(i,j)

        if (dist_ij >0) %there is a delay!
            %HERE, if the sub-systems are not homogeneous, these indeces
            %change and you need to consider how big the other inputs
            %and outputs were
            % (i.e, prev_inputs = sum_{k3=1}^{i} m_k    and
            %       prev_measurements = sum_{k3=1}^{j}p_k     !!!)
            row_start = (i - 1) * m_i + 1;
            row_end = i * m_i;
            col_start = (j - 1) * p_i + 1;
            col_end = j * p_i;
            Q_mask(row_start:row_end,col_start:col_end, 1) = 0; %Q_0 to zero
            if (dist_ij ~=1) %You need to put to zero some values of Q_i with i>0
                % row_start_C = (i - 1) * m_i + 1;
                % row_end_C = i * m_i;
                % col_start_C = (j - 1) * p_i*options.N_tf + 1;
                % col_end_C = j * p_i*options.N_tf;
                % mask = zeros(m_i,p_i*options.N_tf);
                if (dist_ij == Inf) %means that there are some values which can be nonzero
                    % disp("NOT INF!")
                    Q_mask(row_start:row_end,col_start:col_end, 2:end) = 0;
                else
                    Q_mask(row_start:row_end,col_start:col_end, 2:dist_ij) = 0;
                end
            end
        end
    end
end


for k = 1:options.N_tf+1
    [row, col] = find(Q_mask(:,:,k));
    % Create an sdpvar vector with the number of variables equal to the number of non-zero elements in mask
    num_vars = nnz(Q_mask(:,:,k));
    variables = sdpvar(num_vars, 1);
    
    % Create the sparse sdpvar matrix Cq with the sparsity pattern defined by mask
    Q_sdp{k} = sparse(row, col, variables, size(Q_mask(:,:,k), 1), size(Q_mask(:,:,k), 2));

    % Q_sdp(:,:,k) = Q_sdp(:,:,k).*Q_mask(:,:,k);
end


%% LMI Constraints
constraints = [];

%% Objective
n_z_vectors = options.number_points;
random_phase_shift = 0;%rand * 2 * pi;  % Random phase shift in [0, 2*pi]
omega_vector = (0:n_z_vectors-1)/n_z_vectors;  % Evenly spaced angles
z_vector = exp(1i * (2 * pi * omega_vector + random_phase_shift));  % Apply phase shift

i_want_to_grid_more_the_edges = false;
if i_want_to_grid_more_the_edges
    % Add more points between omega = [0, 0.1] and [0.9, 1]
    additional_points = 100; % Number of additional points in each interval
    omega_additional_1 = linspace(0, 0.1, additional_points);
    omega_additional_2 = linspace(0.9, 1, additional_points);
    
    % % Remove the endpoints to avoid duplicates
    % omega_additional_1 = omega_additional_1(2:end-1);
    % omega_additional_2 = omega_additional_2(2:end-1);
    
    % Combine the original and additional points
    omega_vector = [omega_vector, omega_additional_1, omega_additional_2];
    omega_vector = unique(omega_vector);
    z_vector = exp(1i * (2 * pi * omega_vector + random_phase_shift));
    
    % Update the number of z vectors
    n_z_vectors = length(z_vector);
end

%% Adding weights to the objective function on different frequencies
% Parameters
Omega0 = 0.64;       % Resonance at pi/4 (normalized frequency)
r = 0.95;            % Resonance sharpness (adjust as needed)

[A_weight_i, B_weight_i, C_weight_i, D_weight_i] = peak_filter_generator(sys.plant.Ts, Omega0, r, false);
transfer_function_i = ss(A_weight_i, B_weight_i, C_weight_i, D_weight_i, sys.plant.Ts);
transfer_function_i = transfer_function_i / norm(transfer_function_i,2);
weighting_tf = tf(1,1) * eye(sys.n_w); % Identity transfer function for weighting
chosen_channel = 5;
weighting_tf(chosen_channel, chosen_channel) = transfer_function_i; 
disp("DEBUG: WEIGHTING TF AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")


objective = 0 ;
for i=1:n_z_vectors
    z_i = z_vector(i);

    P11_t = evalfr(P11,z_i);
    P12_t = evalfr(P12,z_i);
    P21_t = evalfr(P21,z_i);

    %Trimming matrices
    P11_t = trim_matrix(P11_t,options.trim_tol);
    P12_t = trim_matrix(P12_t,options.trim_tol);
    P21_t = trim_matrix(P21_t,options.trim_tol);
    
    Qsdp_t = evalfr_for_sdpvar(Q_sdp,z_i);
    
    Tzw_tf_t = P11_t - P12_t*Qsdp_t*P21_t;

    weighting_tf_t  = evalfr(weighting_tf, z_i); %DEBUG: TO REMOVE
    weighting_tf_t = trim_matrix(weighting_tf_t, options.trim_tol); %DEBUG: TO REMOVE

    Tzw_tf_t = Tzw_tf_t * weighting_tf_t;

    objective = objective + norm(Tzw_tf_t,'fro');
end

%% Optimizing ...
sol = optimize(constraints, objective, options);
if ~(sol.problem == 0)
    error("Error during calculus of the distributed. Problem status: " +num2str(sol.problem) + " ( " + yalmiperror(sol.problem) + ")");
end
disp("    Problem Solved... Retrieving solutions");
%% Result
Q_sdp_value = zeros(sys.m,sys.p,options.N_tf+1);

for k = 1:options.N_tf+1
    Q_sdp_value(:,:,k) = value(Q_sdp{k});
end

Q_opt = FIR2ss(Q_sdp_value, sys.plant.Ts);
objective = value(objective);
K_opt = get_K_given_Q(sys,Q_opt);
% Q_opt= FIR2tf(value(Q_sdp), sys.plant.Ts); 
end




function [Q_in_z_i] = evalfr_for_sdpvar(sdp_tf,z_i)
%This function evaluates the transfer function at a given point z_i
%It is used to evaluate the transfer function at the point z_i
%when the transfer function is represented as an sdpvar of dimension m x p x (N_tf + 1)

%Input:
%sdp_tf: sdpvar of dimension m x p x (N_tf + 1)
%z_i: complex number

%Output:
%Q_in_z_i: value of the transfer function at the point z_i

N_tf_plus_1 = size(sdp_tf,1);
Q_in_z_i = 1*sdp_tf{1};

for k = 2:N_tf_plus_1
    Q_in_z_i = Q_in_z_i + sdp_tf{k} * z_i^(-k+1);
end
end


