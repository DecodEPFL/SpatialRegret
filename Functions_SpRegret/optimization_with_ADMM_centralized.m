function [K_opt, Phi_solution, objective, miscellanea] = optimization_with_ADMM_centralized(sys, delays_matrix, options)
% OPTIMIZATION_WITH_ADMM_CENTRALIZED
%   ADMM‐based solver for distributed SLS
%
% Inputs:
%   sys           – system struct (n,m,p,C1,D11,D12,D21,B1,…)
%   delays_matrix – network delay pattern
%   options       – struct with fields
%                   .N_tf, .N_FIR_closed_loop, .rho, .tol_ADMM,
%                   .max_ADMM_iter, [.warmstart,.warmStartSolution],…
%
% Outputs:
%   K_opt         – resulting controller (not yet implemented)
%   Phis_spreg    – struct with .mat.R, .M, .N, .L
%   objective     – final ∞-norm cost
%   miscellanea   – residual & cost history

%% 1) Dimensions & horizons
n = sys.n;  m = sys.m;  p = sys.p;
N_tf = options.N_tf;
N_tf_closed_loop = max(options.N_FIR_closed_loop, N_tf);

%% 2) Initialize ADMM variables
Phi_k   = zeros(n+m, (n+p)*(N_tf_closed_loop+1));
Psi_k   = Phi_k;
Lambda_k= Phi_k;

if isfield(options,'warmstart') && options.warmstart
    sol = options.warmStartSolution;
    Phi_k = [sol.R, sol.N; sol.M, sol.L];
    Psi_k = Phi_k;
    Lambda_k = zeros(size(Phi_k));
end

%% 3) Sparsity
L_sparsity      = sls_sparsity_constraints(sys, delays_matrix, N_tf, N_tf_closed_loop);

%% 4) Preallocate history
maxIter = options.max_ADMM_iter;
residual_vector = zeros(maxIter,1);
cost_L1_vector  = zeros(maxIter,1);

%% 5) ADMM loop
rho_k   = options.rho;
tol_ADMM = options.tol_ADMM;

for iter = 1:maxIter
    fprintf('Iter %d: ', iter);
    tic;
    [Phi_k, cost_L1] = row_update_centr(Psi_k, Lambda_k, rho_k, L_sparsity, sys, N_tf_closed_loop, options);
    Psi_k = column_update_centr(Phi_k, Lambda_k, rho_k, L_sparsity, sys, N_tf_closed_loop, options);
    Lambda_k = Lambda_k + Phi_k - Psi_k;
    elapsed = toc;

    [is_converged, residual_k] = check_convergence(Phi_k, Psi_k, tol_ADMM);
    residual_vector(iter) = residual_k;
    cost_L1_vector(iter)  = cost_L1;
    fprintf('(residual: %e) (L1: %e) (ρ_k: %e) (time: %0.2fs)\n', residual_k, cost_L1, rho_k, elapsed);
    if is_converged
        disp('-> converged');
        break;
    end 

    % Increasing rho if you are close to the solution
    if residual_k < 8e-2 && rho_k < 10 
        rho_k = 10 * rho_k;
        disp("updated rho!")
    end
    % Increase more if you are even closer

end

%% 6) Post‐process & return
Rmat = Phi_k(1:n,    1:n*(N_tf_closed_loop+1));
Mmat = Phi_k(n+1:end,1:n*(N_tf_closed_loop+1));
Nmat = Phi_k(1:n,    n*(N_tf_closed_loop+1)+1:end);
Lmat = Phi_k(n+1:end,n*(N_tf_closed_loop+1)+1:end);

[K_opt, Phi_solution] = sls_postprocessing(Rmat,Mmat,Nmat,Lmat,N_tf_closed_loop,sys.plant.Ts);
Phi_solution.mat.R = Rmat; Phi_solution.mat.M = Mmat; Phi_solution.mat.N = Nmat; Phi_solution.mat.L = Lmat;
Phi_solution.mat.Phi = Phi_k;

objective = cost_L1;
miscellanea.residual_vector = residual_vector(1:iter);
miscellanea.cost_L1_vector  = cost_L1_vector(1:iter);
end


%% Row update function
function [Phi_k_sol, cost_L1] = row_update_centr(Psi_k, Lambda_k, rho_k, L_sparsity, sys, N_tf_closed_loop, options)
clear yalmip;
% Row-wise update of the decision variables
n = sys.n;
p = sys.p;
m = sys.m;

% Define the variables
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

Phi_k = [Rvar, Nvar;
    Mvar, Lvar];

% Define the constraints
constraints = [];
ach_row_constr = sls_achievability_constraints_only_rowwise(vars, sys, N_tf_closed_loop);
constraints = [constraints, ach_row_constr];

% Define the objective function

objective_dual = rho_k/2 * norm(Phi_k - Psi_k + Lambda_k , 'fro')^2;

% Adding the L1 objective
Tzw_0 = sys.C1*Rvar(:, 1:n) * sys.B1 + ...
    sys.D12*Mvar(:, 1:n) * sys.B1 + ...
    sys.C1*Nvar(:, 1:p) * sys.D21 + ...
    sys.D12*Lvar(:, 1:p) * sys.D21 + sys.D11;

M_matrix = Tzw_0;

for t= 1:N_tf_closed_loop
    Tzw_t = sys.C1  * Rvar(:, t*n+1:(t+1)*n) * sys.B1 + ...
        sys.D12 * Mvar(:, t*n+1:(t+1)*n) * sys.B1 + ...
        sys.C1  * Nvar(:, t*p+1:(t+1)*p) * sys.D21 + ...
        sys.D12 * Lvar(:, t*p+1:(t+1)*p) * sys.D21;
    M_matrix = [M_matrix, Tzw_t];
end

cost_L1 = norm(M_matrix, 'inf');

objective = objective_dual + cost_L1;

% Solve the optimization problem
% options = sdpsettings('verbose', 0, 'solver', 'gurobi');
diagnostics = optimize(constraints, objective, options);
if diagnostics.problem ~= 0
    error('Row update optimization problem is infeasible or unbounded.');
end
% Extract the solution
Phi_k_sol = value(Phi_k);
cost_L1   = value(cost_L1);
% dual_cost  = value(objective_dual);

end

%% Column update function
function Psi_k_sol = column_update_centr(Phi_k, Lambda_k, rho_k, L_sparsity, sys, N_tf_closed_loop, options)
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

objective = norm(Psi_k - Phi_k - Lambda_k , 'fro')^2;

% Solve the optimization problem
diagnostics = optimize(constraints, objective, options);
if diagnostics.problem ~= 0
    error("Error during Column update!!! Problem status: " +num2str(diagnostics.problem) + " ( " + yalmiperror(diagnostics.problem) + ")");
end
% Extract the solution
Psi_k_sol = value(Psi_k);
end

%% Check convergence function
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