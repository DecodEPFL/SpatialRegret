function [K_opt, Phis_spreg, objective, miscellanea]  =          distributed_optimization_L1_SLS(sys,delays_matrix,options)
clear yalmip;

%% Problem dimensions and horizons
n = sys.n;
m = sys.m;
p = sys.p;
N_tf = options.N_tf;
N_tf_closed_loop = max(options.N_FIR_closed_loop, N_tf);
options.N_tf_closed_loop = N_tf_closed_loop;

%% Preallocate ADMM variables
Phi_k    = zeros(n+m, (n+p)*(N_tf_closed_loop+1));  % row variable
Psi_k    = Phi_k;                                   % column variable
Lambda_k = Phi_k;                                   % dual variable

%% Warm‐start (optional)
if isfield(options, 'warmstart') && options.warmstart
    % options.warmStartSolution has fields R,N; M,L
    sol = options.warmStartSolution;
    initial = [sol.R, sol.N; sol.M, sol.L];
    Phi_k    = initial;
    Psi_k    = initial;
    Lambda_k = Phi_k - Psi_k;
end

%% Sparsity mask for both Phi_k and Psi_k
L_sparsity = sls_sparsity_constraints(sys, delays_matrix, N_tf, N_tf_closed_loop);
mask = [ones(n, n*(N_tf_closed_loop+1)), ones(n, p*(N_tf_closed_loop+1));
        ones(m, n*(N_tf_closed_loop+1)), L_sparsity];
options.Sparsity_mask = mask;

%% ADMM settings
max_ADMM_iter   = options.max_ADMM_iter;
rho_k           = options.rho;
tol_ADMM        = options.tol_ADMM;
row_index_vector = partitionIndices(n+m, options.N_row_indices, 1);

residual_vector = zeros(max_ADMM_iter,1);
cost_L1_vector  = zeros(max_ADMM_iter,1);

%% Main ADMM loop
for iter = 1:max_ADMM_iter
    fprintf('Iter %d: ', iter);
    tic
    [Phi_k, cost_L1] = row_update(Psi_k, Lambda_k, rho_k, sys, row_index_vector, options);
    Psi_k            = column_update_NOT_DISTRIBUTED(Phi_k, Lambda_k, rho_k, L_sparsity, sys, N_tf_closed_loop, options); % Easy projection, not distributed (can be made distributed if needed)
    Lambda_k         = Lambda_k + Phi_k - Psi_k;
    elapsed = toc;

    [is_converged, residual_k] = check_convergence(Phi_k, Psi_k, tol_ADMM);
    residual_vector(iter) = residual_k;
    cost_L1_vector(iter)  = cost_L1;
    fprintf('(residual: %e) (L1 cost: %e) (rho_k: %e) (time: %0.2fs)\n', residual_k, cost_L1, rho_k, elapsed);
    if is_converged
        fprintf('Converged at iteration %d\n', iter);
        break;
    end  
end

%% Collect outputs
K_opt = [];  % DEBUG: controller synthesis not yet implemented
Phis_spreg.mat.R = Phi_k(1:n,                1:n*(N_tf_closed_loop+1));
Phis_spreg.mat.M = Phi_k(n+1:end,            1:n*(N_tf_closed_loop+1));
Phis_spreg.mat.N = Phi_k(1:n,                n*(N_tf_closed_loop+1)+1:end);
Phis_spreg.mat.L = Phi_k(n+1:end,            n*(N_tf_closed_loop+1)+1:end);

objective    = cost_L1;
miscellanea.residual_vector = residual_vector(1:iter);
miscellanea.cost_L1_vector  = cost_L1_vector(1:iter);
end




function [Phi_k_sol, cost_L1] = row_update(Psi_k, Lambda_k, rho_k, sys, row_index_vector, options)
% Row update of the decision variables. It uses the golden ratio search to find the optimal gamma.

Phi_k_sol = zeros(size(Psi_k)); % Initialize the output variable

% Extract the hyperparameters
max_iter_gr = options.max_iter_gr;
gamma_a = options.gamma_lower;
gamma_b = options.gamma_upper;
tol_gr = options.tol_gr;

% Golden ratio search
r = (3-sqrt(5))/2; % golden ratio

gamma_c = gamma_a + r*(gamma_b-gamma_a);        % preparation
gamma_d = gamma_b - r*(gamma_b-gamma_a);

clear yalmip;
size_of_Phi_k = size(Psi_k);


row_optimizers = row_optimizer_generator(Psi_k - Lambda_k, rho_k, sys, options.N_tf_closed_loop, row_index_vector, options.Sparsity_mask, options);


[fc, Phi_c] = eval_row_subproblems(gamma_c, row_optimizers, row_index_vector, size_of_Phi_k);
[fd, Phi_d] = eval_row_subproblems(gamma_d, row_optimizers, row_index_vector, size_of_Phi_k);

assert(fd < Inf, 'gamma_upper is infeasible. Choose a larger value for gamma_upper.');

iter = 0;

while abs(gamma_b - gamma_a) >= tol_gr && iter < max_iter_gr
    iter = iter + 1;
    % Print current interval and test points
    % fprintf('    GR %2d: [a,b]=[%.3f,%.3f], c=%.3f(f=%.3f), d=%.3f(f=%.3f)\n', iter, gamma_a, gamma_b, gamma_c, fc, gamma_d, fd);

    if fc < fd
        % Store the solution
        Phi_k_sol = Phi_c;
        cost_L1 = gamma_c;

        gamma_b = gamma_d;
        % shift d into c
        gamma_d = gamma_c;
        fd = fc;
        % new c
        gamma_c = gamma_a + r * (gamma_b - gamma_a);
        % fprintf('     new c = %e\n', gamma_c);
        [fc, Phi_c] = eval_row_subproblems(gamma_c, row_optimizers, row_index_vector, size_of_Phi_k);
        % fprintf('     evaluated fc at c = %e -> fc = %e\n', gamma_c, fc);
    else
        % Store the solution
        Phi_k_sol = Phi_d;
        cost_L1 = gamma_d;

        gamma_a = gamma_c;
        % shift c into d
        gamma_c = gamma_d;
        fc = fd;
        % new d
        gamma_d = gamma_b - r * (gamma_b - gamma_a);
        % fprintf('     new d = %e\n', gamma_d);
        [fd, Phi_d] = eval_row_subproblems(gamma_d, row_optimizers, row_index_vector, size_of_Phi_k);
        % fprintf('     evaluated fd at d = %e -> fd = %e\n', gamma_d, fd);
    end

end
% fprintf('    Distributed method done in %d iters --> gamma=%.4f\n', iter, cost_L1);
end





function [total_objective, Phi_k_temp] = eval_row_subproblems(gamma_i, row_optimizers, row_index_vector, size_of_Phi)
numRows = numel(row_optimizers);
Phi_k_temp = zeros(size_of_Phi);  % Initialize the output variable
total_dual = 0; % Initialize the total dual variable
% Loop through each row subproblem
for i = 1:numRows
    idx = row_index_vector(i,1):row_index_vector(i,2);

    % call the optimizer
    [solution_i, problem_flag_i] = row_optimizers{i}(gamma_i);  % solution_i  = cell { dual_i , Phi_k_row }

    % Check for infeasibility
    if problem_flag_i ~= 0
        total_objective = Inf; return % Infeasible so you don't need to check the other rows
    end
    dual_i = solution_i{1};
    total_dual = total_dual + dual_i;
    Phi_k_temp(idx,:) = solution_i{2};
end
total_objective = total_dual + gamma_i;
end



function Psi_k_sol = column_update_NOT_DISTRIBUTED(Phi_k, Lambda_k, rho, L_sparsity, sys, N_tf_closed_loop, options)
clear yalmip;
% Row-wise update of the decision variables
n = sys.n;
p = sys.p;
m = sys.m;

Rvar = sdpvar(n,n*(N_tf_closed_loop+1),'full');          % decision variables for R
Mvar = sdpvar(m,n*(N_tf_closed_loop+1),'full');          % decision variables for M
Nvar = sdpvar(n,p*(N_tf_closed_loop+1),'full');          % decision variables for N

%Sparsity Constraints
Lvar = sdpvar_sparse(L_sparsity); % Lvar = sdpvar(m,p*(N_tf+1));          % decision variables for L

vars = struct;
vars.Rvar = Rvar;
vars.Mvar = Mvar;
vars.Nvar = Nvar;
vars.Lvar = Lvar;

Psi_k = [Rvar, Nvar;
    Mvar, Lvar];

% Define the constraints
constraints = [];
ach_col_constr = sls_achievability_constraints_only_columnwise(vars, sys, N_tf_closed_loop);
constraints = [constraints, ach_col_constr];

% Define the objective function    Objective = rho/2*||Psi_k - Phi_k - Lambda_k||^2_{H2}
% objective = rho/2 * norm(Psi_k - Phi_k - Lambda_k , 'fro')^2;
objective = norm(Psi_k - Phi_k - Lambda_k , 'fro')^2;

% Solve the optimization problem
% options = sdpsettings('verbose', 0, 'solver', 'mosek');
diagnostics = optimize(constraints, objective, options);
if diagnostics.problem ~= 0
    error("Error during Column update!!! Problem status: " +num2str(diagnostics.problem) + " ( " + yalmiperror(diagnostics.problem) + ")");
end
% Extract the solution
Psi_k_sol = value(Psi_k);
end



function index_vector = partitionIndices(total_elements, num_groups, dimension)
if nargin < 3
    dimension = 1;
end
assert(num_groups<=total_elements, 'Number of groups must be less than or equal to total elements.');
% Compute the base number of elements per group.
base_size = floor(total_elements / num_groups);
% Preallocate the result array.
index_vector = zeros(num_groups, 2);

% Set the starting index.
current_index = 1;

for i = 1:num_groups
    if i < num_groups
        % For groups 1 to num_groups-1, assign exactly base_size elements.
        group_size = base_size;
    else
        % The last group gets all the remaining elements.
        group_size = total_elements - current_index + 1;
    end
    % Store first and last indices for this group.
    index_vector(i, :) = [current_index, current_index + group_size - 1];
    % Update the current index.
    current_index = current_index + group_size;
end

% Adjust the output dimension.
if dimension == 2
    index_vector = index_vector';
end
end


function [converged, residual] = check_convergence(Phi_k, Psi_k, tol)
% Check convergence criteria for ADMM
% This function checks if the primal residual is below the specified tolerance.
% Inputs:
%   Phi_k: Current value of the primal variable
%   Psi_k: Current value of the dual variable
%   tol: Tolerance for convergence
% Outputs:
%   converged: Boolean indicating whether the algorithm has converged

% Compute ADMM primal residual
residual = norm(Phi_k - Psi_k, 'fro');
% fprintf('(residual = %e)', r);
converged = (residual <= tol);
end


function cost_L1 = compute_cost_L1(Phi_k, sys, N_tf_closed_loop)
% Phi_k = [R  N;  M  L] stacked horizontally in time
n = sys.n;  p = sys.p;
% zero‑lag
Rmat = Phi_k(1:n, 1:n*(N_tf_closed_loop+1));
Mmat = Phi_k(n+1:end, 1:n*(N_tf_closed_loop+1));
Nmat = Phi_k(1:n, n*(N_tf_closed_loop+1)+1:end);
Lmat = Phi_k(n+1:end, n*(N_tf_closed_loop+1)+1:end);

R0 = Rmat(:, 1:n);
M0 = Mmat(:, 1:n);
N0 = Nmat(:, 1:p);
L0 = Lmat(:, 1:p);

T0 = sys.C1*R0*sys.B1 ...
    + sys.D12*M0*sys.B1 ...
    + sys.C1*N0*sys.D21 ...
    + sys.D12*L0*sys.D21 ...
    + sys.D11;
Matrix_cost = T0;

for t = 1:N_tf_closed_loop
    Rt = Rmat(:, t*n+1:(t+1)*n);
    Mt = Mmat(:, t*n+1:(t+1)*n);
    Nt = Nmat(:, t*p+1:(t+1)*p);
    Lt = Lmat(:, t*p+1:(t+1)*p);

    Tt = sys.C1*Rt*sys.B1 ...
        + sys.D12*Mt*sys.B1 ...
        + sys.C1*Nt*sys.D21 ...
        + sys.D12*Lt*sys.D21;
    Matrix_cost = [Matrix_cost, Tt];
end
cost_L1 = norm(Matrix_cost,'inf');
% cost_L1 = max(sum(abs(Matrix_cost),2));
end


function dual_cost = compute_dual_cost(Phi_k, Psi_k, Lambda_k, rho)
% Centralized dual cost:  (ρ/2) * ||Φ - Ψ + Λ||_F^2
diff = Phi_k - Psi_k + Lambda_k;
dual_cost = rho/2 * norm(diff, 'fro')^2;
end


function total_cost = compute_total_cost(Phi_k, Psi_k, Lambda_k, rho, sys, N_tf_closed_loop)
% Compute the total cost as the sum of the L1 cost and the dual cost
L1_cost   = compute_cost_L1  (Phi_k, sys, N_tf_closed_loop);
dual_cost = compute_dual_cost(Phi_k, Psi_k, Lambda_k, rho);
total_cost = L1_cost + dual_cost;
end

function [rho_new, Lambda] = update_rho_rb(Phi, Psi, Psi_prev, rho, Lambda, mu, tau)
%UPDATE_RHO_RB  Residual‐Balancing ADMM rho‐update
%
%   [rho_new, Lambda] = update_rho_rb(Phi, Psi, Psi_prev, rho, Lambda, mu, tau)
%
%   Inputs:
%     Phi      – current Phi^(k)
%     Psi      – current Psi^(k)
%     Psi_prev – previous Psi^(k-1)
%     rho      – previous rho^(k-1)
%     Lambda   – current scaled dual Lambda^(k) (before rescaling)
%     mu       – balance tolerance (e.g. 10)
%     tau      – scaling factor (e.g. 2)
%
%   Outputs:
%     rho_new  – updated penalty parameter
%     Lambda   – rescaled dual variable

% 1) compute residuals
r = Phi - Psi;
s = rho * (Psi - Psi_prev);

nr = norm(r, 'fro');
ns = norm(s, 'fro');

% 2) update rho
if nr > mu * ns
    rho_new = rho * tau;
elseif ns > mu * nr
    rho_new = rho / tau;
else
    rho_new = rho;
end

% 3) rescale Lambda
Lambda = (rho / rho_new) * Lambda;
end

function [rho_new, Lambda] = update_rho_spec(r, r_prev, Lambda, Lambda_prev, rho, rho_min, rho_max)
%UPDATE_RHO_SPEC  Spectral (Barzilai–Borwein) ADMM rho‐update
%
%   [rho_new, Lambda] = update_rho_spec(r, r_prev, Lambda, Lambda_prev, rho, rho_min, rho_max)
%
%   Inputs:
%     r            – current primal residual (Phi - Psi)
%     r_prev       – previous primal residual
%     Lambda       – current scaled dual Lambda^(k)
%     Lambda_prev  – previous scaled dual Lambda^(k-1)
%     rho          – previous rho^(k-1)
%     rho_min      – lower safeguard (e.g. 1e-6)
%     rho_max      – upper safeguard (e.g. 1e+6)
%
%   Outputs:
%     rho_new      – updated penalty parameter
%     Lambda       – rescaled dual variable

% 1) differences
dr      = r       - r_prev;
dLambda = Lambda  - Lambda_prev;

% 2) Barzilai–Borwein estimate
num = dr(:)' * dr(:);
den = dr(:)' * dLambda(:);

if den <= 0
    rho_bar = rho;    % fallback
else
    rho_bar = num / den;
end

% 3) safeguard
rho_new = min( max(rho_bar, rho_min), rho_max );

% 4) rescale Lambda
Lambda = (rho / rho_new) * Lambda;
end

function row_optimizers = row_optimizer_generator(Psi_minus_Lambda_k, rho_k, sys, N_tf_closed_loop, row_index_vector, Sparsity_mask, options)
% Builds one YALMIP optimizer per row‐block.
n = sys.n; p = sys.p;

numRowSubproblems = size(row_index_vector,1);
row_optimizers = cell(numRowSubproblems,1);
% fprintf('\n    Building %d row optimizers... ', numRowSubproblems);
for i = 1:numRowSubproblems
    
    idx = row_index_vector(i,1):row_index_vector(i,2); % Indices for this row
    
    Phi_row = sdpvar_sparse(Sparsity_mask(idx,:)); % decision variable for this row
    
    gamma_par  = sdpvar(1,1);  % parameters

    % build the optimizer
    %
    % achievability constraints for this row
    Constraints = sls_achievability_constraints_distrib(Phi_row, sys, N_tf_closed_loop, 'row', row_index_vector(i,:));

    % dual term
    % obj_dual = rho/2 * norm(Phi_row - Psi_row + Lambda_row, 'fro')^2;
    obj_dual = rho_k/2 * norm(Phi_row - Psi_minus_Lambda_k(idx,:), 'fro')^2;

    % build M matrix for L∞ constraint
    LHS = [sys.C1, sys.D12];
    WeightedPhi_row = LHS(idx,idx)*Phi_row;
    Mmat = WeightedPhi_row(:,1:n)*sys.B1 + WeightedPhi_row(:,n*(N_tf_closed_loop+1)+[1:p])*sys.D21 + sys.D11(idx,:);
    for t=1:N_tf_closed_loop
        WeightedPhi_LEFT = WeightedPhi_row(:,t*n+1:(t+1)*n);
        WeightedPhi_RIGHT = WeightedPhi_row(:,n*(N_tf_closed_loop+1)+(t*p+1:(t+1)*p));
        Mmat = [Mmat, WeightedPhi_LEFT*sys.B1 + WeightedPhi_RIGHT*sys.D21];
    end
    Constraints = [Constraints, norm(Mmat,'inf') <= gamma_par];

    % build optimizer
    row_optimizers{i} = optimizer(Constraints, obj_dual, options, ...
        {gamma_par}, ...      % Parameters (Inputs)
        {obj_dual, Phi_row}); %Outputs
    % fprintf('->')
end
% fprintf('done.');
end