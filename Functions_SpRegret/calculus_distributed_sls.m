function [K_opt, Phis, objective] = calculus_distributed_sls(sys,problem, delays_matrix,options)
n  = sys.n;
p = sys.p;
m = sys.m;
N_tf= options.N_tf;       % FIR horizon

N_tf_closed_loop = max( options.N_FIR_closed_loop , N_tf ); % FIR horizon for closed-loop system. Checking that it is not smaller than N_tf
%% Defining Optimization Variables
Rvar = sdpvar(n,n*(N_tf_closed_loop+1),'full');          % decision variables for R
Mvar = sdpvar(m,n*(N_tf_closed_loop+1),'full');          % decision variables for M
Nvar = sdpvar(n,p*(N_tf_closed_loop+1),'full');          % decision variables for N

%Sparsity Constraints
L_sparsity = sls_sparsity_constraints(sys,delays_matrix,N_tf,N_tf_closed_loop);
Lvar = sdpvar_sparse(L_sparsity); % Lvar = sdpvar(m,p*(N_tf+1));          % decision variables for L

vars = struct;
vars.Rvar = Rvar;
vars.Mvar = Mvar;
vars.Nvar = Nvar;
vars.Lvar = Lvar;

%% Constraints
constraints = [];

%Achievability Constraints
[achievability_constraints] = sls_achievability_constraints(vars, sys, N_tf_closed_loop);
constraints = [constraints;achievability_constraints];

%% Cost
[objective, cost_constraints]  = SLS_cost(vars,sys,options,problem);
constraints = [constraints;cost_constraints];

% Uncomment and modify warm-start assignments if a saved solution exists.
if isfield(options,'warmstart') && options.warmstart
    % options.warmStartSolution is assumed to be a struct with fields L, M, N, R
    warmstart(Lvar, options.warmStartSolution.L);
    warmstart(Mvar, options.warmStartSolution.M);
    warmstart(Nvar, options.warmStartSolution.N);
    warmstart(Rvar, options.warmStartSolution.R);
end

%% Solving the Problem
sol = optimize(constraints, objective, options);

if ~(sol.problem == 0)
    error("Error during synthesis of the distributed controller! Problem status: " +num2str(sol.problem) + " ( " + yalmiperror(sol.problem) + ")");
end
disp("    Problem Solved... Retrieving solutions");

%Storing results...
Rmat = value(Rvar); Mmat = value(Mvar); Nmat = value(Nvar); Lmat = value(Lvar);
% Rmat = trim_matrix(Rmat); Mmat = trim_matrix(Mmat); Nmat = trim_matrix(Nmat); Lmat = trim_matrix(Lmat); % Trimming numerical errors from the matrices

%  closed-loop responses from optimization

[K_opt, Phis] = sls_postprocessing(Rmat,Mmat,Nmat,Lmat,N_tf_closed_loop,sys.plant.Ts);
Phis.mat.R = Rmat; Phis.mat.M = Mmat; Phis.mat.N = Nmat; Phis.mat.L = Lmat;
Phis.mat.Phi = [Rmat Nmat; Mmat Lmat];

if strcmp(problem,'h2')
    objective = sqrt(value(objective));
else
    objective = value(objective);
end