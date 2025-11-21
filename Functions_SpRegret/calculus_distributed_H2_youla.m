function [Q_opt, H2_norm] = calculus_distributed_H2_youla(sys, delays_matrix,options)

P11 = sys.P11;
P12 = sys.P12;
P21 = sys.P21;
%Building static part of the Youla Parameter Q (i.e., Aq Bq)
Aq_i = diag(ones(options.N_tf-1,1),1); % The Aq matrix for each agent!
Aq = kron(eye(sys.p), Aq_i);

Bq_i = zeros(options.N_tf,1); Bq_i(end,1) = 1; %Bq_i = [0 ... 0 1]'
Bq = kron(eye(sys.p), Bq_i);

%Computing the HAT_matrices
A1_hat = blkdiag(P11.A,P12.A);
A2_hat = [P21.A, zeros(size(P21.A,1), size(Aq,2));
    Bq*P21.C, Aq];

B11_hat = [P11.B; zeros(size(P12.B,1), size(P11.B,2))];
B_hat = [zeros(size(P11.B,1), size(P12.B,2)); P12.B];
B21_hat = [P21.B; Bq*P21.D];

C11_hat = [P11.C, -P12.C];
C_hat = [zeros(sys.p*options.N_tf, size(P21.C,2)), eye(sys.p*options.N_tf);
    P21.C, zeros(size(P21.C,1), sys.p*options.N_tf)];

D11_hat = P11.D;
D12_hat = -P12.D;
D21_hat = [zeros(sys.p*options.N_tf,size(P11.D,2));P21.D];

%% Optimization Variables
gam = sdpvar(1);
E = sdpvar(size(A1_hat,1),size(A1_hat,1), 'symmetric');
R = sdpvar(size(A2_hat,1),size(A2_hat,1), 'symmetric');
S = sdpvar(size(E,1),size(R,2),'full');
% Cq = sdpvar(sys.m,sys.p*options.N_tf,'full');
% Dq = sdpvar(sys.m,sys.p,'full');
Z = sdpvar(size(D11_hat,1),size(D11_hat,1),'symmetric');

%% Constraints
constraints = [];

%Information Sparsities
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


%Bar variables
Q_static = [Cq Dq];
X_bar = blkdiag(E,R);
A_bar = [A1_hat*E,       A1_hat*S+B_hat*Q_static*C_hat - S*A2_hat;
    zeros(size(R,1),size(E,2)), R*A2_hat];
B_bar = [B11_hat+B_hat*Q_static*D21_hat-S*B21_hat;
    R*B21_hat];
C_bar = [C11_hat*E, C11_hat*S+D12_hat*Q_static*C_hat];
D_bar = D11_hat + D12_hat*Q_static*D21_hat;

%LMIs
Matrix1 = [X_bar, zeros(size(X_bar,1), size(B_bar,2)), A_bar';
    zeros(size(B_bar,2),size(X_bar,1)), gam*eye(size(B_bar,2)), B_bar';
    A_bar, B_bar, X_bar];

Matrix2 = [X_bar, zeros(size(X_bar,1), size(D_bar,2)), C_bar';
    zeros(size(D_bar,2),size(X_bar,1)), gam*eye(size(D_bar,2)), D_bar';
    C_bar, D_bar, Z];

constraints = [constraints; Matrix1 >= 0; Matrix2 >= 0; trace(Z)<= gam]; 

% t = sdpvar(1);disp("Added some slack variables for debugging!!")
% LMI_total = blkdiag(Matrix1,Matrix2);
% constraints = [constraints; LMI_total >=t*eye(size(LMI_total))  trace(Z)<= gam]; disp("Modified the constraints");

%% Objective
objective = gam;
% objective = 0; disp("Changed the objective cost for debugging!");

%% Optimizing ...
sol = optimize(constraints, objective, options);

if sol.problem == 4
    warning("The problem didn't end properly. Problem status 4 was returned. The solutions may not be feasible.")
elseif sol.problem == 0
    disp("    Problem Solved... Retrieving solutions"); 
else
    error("Error during calculus of the distributed. Problem status: " +num2str(sol.problem) + " ( " + yalmiperror(sol.problem) + ")");
end

%% Result
Cq = value(Cq);
Dq = value(Dq);
Cq = trim_matrix(Cq); Dq = trim_matrix(Dq);
Q_opt = ss(Aq,Bq,Cq,Dq,sys.plant.Ts);
H2_norm = value(gam);
end

