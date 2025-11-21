function [K_opt, Phis_or_Q, objective] = spatial_regret_weighted(sys,oracle_LFT,delays_matrix,options)
clear yalmip;
method = lower(options.method);
disp("SPATIAL REGRET WEIGHTEEED");

if strcmp(method,'sls')
    [K_opt, Phis_or_Q, objective] = SLS_spatial_regret_sampled(sys,oracle_LFT,delays_matrix,options);

elseif strcmp(method,'sampled_youla')
    [K_opt, Phis_or_Q, objective] = Youla_spatial_regret_sampled(sys,oracle_LFT,delays_matrix,options);

elseif strcmp(method,'tf_sampled')
    [K_opt, Phis_or_Q, objective] = Tf_spatial_regret_sampled(sys,oracle_LFT,delays_matrix,options);

else
    error("Method not recognized. Please choose between 'SLS' and 'sampled_youla' or 'tf_sampled'.");
end

end



function [K_opt, Phis, objective] = SLS_spatial_regret_sampled(sys,oracle_LFT,delays_matrix,options)

n  = sys.n;
p = sys.p;
m = sys.m;
N_tf     = options.N_tf;       % FIR horizon

%% Defining Optimization Variables
Rvar = sdpvar(n,n*(N_tf+1));          % decision variables for R
Mvar = sdpvar(m,n*(N_tf+1));          % decision variables for M
Nvar = sdpvar(n,p*(N_tf+1));          % decision variables for N
Lvar = sdpvar(m,p*(N_tf+1));          % decision variables for L

%% Constraints
constraints = [];

%Achievability Constraints
achievability_constraints = sls_achievability_constraints(Rvar,Mvar,Nvar,Lvar,sys,N_tf);
constraints = [constraints;achievability_constraints];

%Sparsity Constraints
sparsity_constraints = sls_sparsity_constraints(sys,Lvar,delays_matrix,N_tf);
constraints = [constraints;sparsity_constraints];

% Cost
lambda = sdpvar(1);
objective = lambda;
constraints = [constraints, lambda>=0];
n_z_vectors = options.number_points;
z_vector = exp(1i*2*pi*(0:n_z_vectors-1)/n_z_vectors);

for i=1:n_z_vectors
    z_i = z_vector(i);
    Tzw_0 = sys.D12*Lvar(:,1:p)*sys.D21 + sys.D11;

    Tzw_tf_t = Tzw_0;
    Oracle_tf_t = evalfr(oracle_LFT,z_i);

    for t = 1:N_tf
        Tzw_t = sys.C1*Rvar(:,t*n+1:(t+1)*n)*sys.B1 + sys.D12*Mvar(:,t*n+1:(t+1)*n)*sys.B1 + sys.C1*Nvar(:,t*p+1:(t+1)*p)*sys.D21 + sys.D12*Lvar(:,(t)*p+1:(t+1)*p)*sys.D21;
        Tzw_tf_t = Tzw_tf_t + z_i^(-t)*Tzw_t;
    end
    LMI_t = [1*eye(sys.n_z), Tzw_tf_t; Tzw_tf_t', lambda*eye(sys.n_w) + Oracle_tf_t'*Oracle_tf_t];
    constraints = [constraints, LMI_t >= 0];
end

% Get a solution
sol = optimize(constraints, objective, options);
if ~(sol.problem == 0)
    error("Error during synthesis of the distributed controller! Problem status: " +num2str(sol.problem) + " ( " + yalmiperror(sol.problem) + ")");
end
disp("    Problem Solved... Retrieving solutions");

%Storing results...
Rmat = value(Rvar); Mmat = value(Mvar); Nmat = value(Nvar); Lmat = value(Lvar);
% Rmat = trim_matrix(Rmat); Mmat = trim_matrix(Mmat); Nmat = trim_matrix(Nmat); Lmat = trim_matrix(Lmat); % Trimming numerical errors from the matrices

%  closed-loop responses from optimization

[K_opt, Phis] = sls_postprocessing(Rmat,Mmat,Nmat,Lmat,N_tf,sys.plant.Ts);
Phis.mat.R = Rmat; Phis.mat.M = Mmat; Phis.mat.N = Nmat; Phis.mat.L = Lmat;

objective = value(objective);
end



function [K_opt, Q_opt, objective] = Youla_spatial_regret_sampled(sys,oracle_LFT,delays_matrix,options)

P11 = sys.P11;
P12 = sys.P12;
P21 = sys.P21;
%Building static part of the Youla Parameter Q (i.e., Aq Bq)
Aq_i = diag(ones(options.N_tf-1,1),1); % The Aq matrix for each agent!
Aq = kron(eye(sys.p), Aq_i);

Bq_i = zeros(options.N_tf,1); Bq_i(end,1) = 1; %Bq_i = [0 ... 0 1]'
Bq = kron(eye(sys.p), Bq_i);




%% Optimization Variables
lambda = sdpvar(1);
% Cq = sdpvar(sys.m,sys.p*options.N_tf,'full'); %Defined Cq later as as sparse variable
% Dq = sdpvar(sys.m,sys.p,'full');              %Defined Dq later as as sparse variable

%% Information Sparsities
N_agents = size(delays_matrix,1);
m_i = sys.m/N_agents;
p_i = sys.p/N_agents;
Cq_mask = ones(sys.m,sys.p*options.N_tf);
Dq_mask = ones(sys.m,sys.p);

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
            Dq_mask(row_start:row_end,col_start:col_end) = 0; %D to zero
            if (dist_ij ~=1) %You need to put to zero some values of the Cq
                row_start_C = (i - 1) * m_i + 1;
                row_end_C = i * m_i;
                col_start_C = (j - 1) * p_i*options.N_tf + 1;
                col_end_C = j * p_i*options.N_tf;
                mask = zeros(m_i,p_i*options.N_tf);
                if (dist_ij ~= Inf) %means that there are some values which can be nonzero
                    % disp("NOT INF!")
                    mask(1:m_i,1:p_i*(options.N_tf-dist_ij+1)) =1;
                end
                % fprintf("Constraint over Cq(%d : %d,%d : %d) == 0 \n\n",row_start_C,row_end_C,col_start_C,col_end_C);
                Cq_mask(row_start_C:row_end_C,col_start_C:col_end_C)= mask;
            end
        end
    end
end
Cq = sdpvar_sparse(Cq_mask);
Dq = sdpvar_sparse(Dq_mask);


%% LMI Constraints
constraints = lambda>=0;
n_z_vectors = options.number_points;
random_phase_shift = 0;%rand * 2 * pi;  % Random phase shift in [0, 2*pi]
omega_vector = (0:n_z_vectors-1)/n_z_vectors;  % Evenly spaced angles
z_vector = exp(1i * (2 * pi * omega_vector + random_phase_shift));  % Apply phase shift


% disp("Adding the frequency peak manually!!!");
% [~, fpeak] = getPeakGain(sys.plant(1,size(delays_matrix,1)));
% z_vector = [z_vector, exp(1i*fpeak*sys.plant.Ts)];
% n_z_vectors = n_z_vectors+1;

for i=1:n_z_vectors
    z_i = z_vector(i);

    Oracle_tf_t = evalfr(oracle_LFT,z_i);

    P11_t = evalfr(P11,z_i);
    P12_t = evalfr(P12,z_i);
    P21_t = evalfr(P21,z_i);

    Tzw_tf_t = P11_t - P12_t*Cq/(z_i*eye(size(Aq,1)) - Aq)*Bq*P21_t -P12_t*Dq*P21_t;

    LMI_t = [1*eye(sys.n_z), Tzw_tf_t; Tzw_tf_t', lambda*eye(sys.n_w) + Oracle_tf_t'*Oracle_tf_t];
    constraints = [constraints, LMI_t >= 0];
end

%% Objective
objective = lambda;

%% Optimizing ...
sol = optimize(constraints, objective, options);
if ~(sol.problem == 0)
    error("Error during calculus of the distributed. Problem status: " +num2str(sol.problem) + " ( " + yalmiperror(sol.problem) + ")");
end
disp("    Problem Solved... Retrieving solutions");
%% Result
Cq = value(Cq);
Dq = value(Dq);
Cq = trim_matrix(Cq); Dq = trim_matrix(Dq);
Q_opt = ss(Aq,Bq,Cq,Dq,sys.plant.Ts);
objective = value(objective);
K_opt = get_K_given_Q(sys,Q_opt); %TODO: Fix this function to work properly... K is not being computed correctly
end



function [K_opt, Q_opt, objective] = Tf_spatial_regret_sampled(sys,oracle_LFT,delays_matrix,options)
P11 = sys.P11;
P12 = sys.P12;
P21 = sys.P21;

%% Optimization Variables
lambda = sdpvar(1);
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
end

%% LMI Constraints
constraints = lambda>=0;

n_z_vectors = options.number_points;
random_phase_shift = 0;%rand * 2 * pi;  % Random phase shift in [0, 2*pi]
omega_vector = (0:n_z_vectors-1)/n_z_vectors;  % Evenly spaced angles
z_vector = exp(1i * (2 * pi * omega_vector + random_phase_shift));  % Apply phase shift

weighted_frequencies = options.weighted_frequencies;
for i=1:n_z_vectors
    z_i = z_vector(i);

    P11_t = evalfr(P11,z_i);
    P12_t = evalfr(P12,z_i);
    P21_t = evalfr(P21,z_i);
    Oracle_tf_t  = evalfr(oracle_LFT,z_i);

    P11_t = trim_matrix(P11_t,options.trim_tol);
    P12_t = trim_matrix(P12_t,options.trim_tol);
    P21_t = trim_matrix(P21_t,options.trim_tol);

    Qsdp_t = evalfr_for_sdpvar(Q_sdp,z_i);

    Tzw_tf_t = P11_t - P12_t*Qsdp_t*P21_t;

    % Finding the frequencies between [0.605,0.83] and [0.605 + pi,0.83 +pi]
    omega_i = omega_vector(i)*2*pi;
    if (omega_i>=weighted_frequencies(1) && omega_i<=weighted_frequencies(2)) || ...
            (omega_i>=weighted_frequencies(1)+pi && weighted_frequencies(1)<=0.83+pi)
        alpha = 0.75;
    else
        alpha = 1;
    end

    LMI_t = [1*eye(sys.n_z), Tzw_tf_t; Tzw_tf_t', lambda*eye(sys.n_w) + alpha*Oracle_tf_t'*Oracle_tf_t];
    constraints = [constraints, LMI_t >= 0];
end

%% Objective
objective = lambda;


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