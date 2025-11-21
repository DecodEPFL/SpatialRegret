function [K_opt, Q_opt, objective] = calculus_distributed_H2_youla_sampled(sys, delays_matrix,options)

P11 = sys.P11;
P12 = sys.P12;
P21 = sys.P21;
%Building static part of the Youla Parameter Q (i.e., Aq Bq)
Aq_i = diag(ones(options.N_tf-1,1),1); % The Aq matrix for each agent!
Aq = kron(eye(sys.p), Aq_i);

Bq_i = zeros(options.N_tf,1); Bq_i(end,1) = 1; %Bq_i = [0 ... 0 1]'
Bq = kron(eye(sys.p), Bq_i);




%% Optimization Variables
% lambda = sdpvar(1);
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
constraints = [];

%% Objective
n_z_vectors = options.number_points;
random_phase_shift = 0;  %rand * 2 * pi;  % Random phase shift in [0, 2*pi]
omega_vector = (0:n_z_vectors-1)/n_z_vectors;  % Evenly spaced angles
z_vector = exp(1i * (2 * pi * omega_vector + random_phase_shift));  % Apply phase shift

objective = 0 ;
eye_q = eye(size(Aq,1));
for i=1:n_z_vectors
    z_i = z_vector(i);

    P11_t = evalfr(P11,z_i);
    P12_t = evalfr(P12,z_i);
    P21_t = evalfr(P21,z_i);
    
    P11_t = trim_matrix(P11_t,options.trim_tol);
    P12_t = trim_matrix(P12_t,options.trim_tol);
    P21_t = trim_matrix(P21_t,options.trim_tol);
    
    Tzw_tf_t = P11_t - P12_t*(Cq/(z_i*eye_q - Aq)*Bq + Dq)*P21_t;
    objective = objective + norm(Tzw_tf_t,'fro')^2;
end

%% Optimizing ...
disp("    Starting Optimization.");
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

