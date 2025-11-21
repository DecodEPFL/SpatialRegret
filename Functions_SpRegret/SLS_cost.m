function [cost, constraints] = SLS_cost(vars, sys, options, problem)
% SLS_cost computes the cost and constraints for the given problem.
% It expects vars to be a structure holding:
%   vars.Rvar, vars.Mvar, vars.Nvar, vars.Lvar
%
% You can also include additional variables in vars if needed.

% Extrapolate variables from the input structure
Rvar = vars.Rvar;
Mvar = vars.Mvar;
Nvar = vars.Nvar;
Lvar = vars.Lvar;


p = sys.p;
n = sys.n;
m = sys.m;
N_tf = max(options.N_tf, options.N_FIR_closed_loop); % FIR horizon
constraints = [];

if strcmp(problem, 'h2')
    cost = 0;
    Tzw_0 = sys.D12 * Lvar(:,1:p) * sys.D21 + sys.D11;
    cost = cost + trace(Tzw_0' * Tzw_0);
    for t = 1:N_tf
        Tzw_t = sys.C1 * Rvar(:, t*n+1:(t+1)*n) * sys.B1 + ...
                sys.D12 * Mvar(:, t*n+1:(t+1)*n) * sys.B1 + ...
                sys.C1 * Nvar(:, t*p+1:(t+1)*p) * sys.D21 + ...
                sys.D12 * Lvar(:, t*p+1:(t+1)*p) * sys.D21;
        cost = cost + trace(Tzw_t' * Tzw_t);
    end
elseif strcmp(problem, 'hinf')
    disp("New Hinf Method")
    lambda = sdpvar(1);
    cost = lambda;
    constraints = [constraints, lambda >= 0];
    n_z_vectors = options.number_points;
    z_vector = exp(1i*2*pi*(0:n_z_vectors-1)/n_z_vectors);

    for i = 1:n_z_vectors
        z_i = z_vector(i);
        Tzw_0 = sys.D12 * Lvar(:,1:p) * sys.D21 + sys.D11;
        Tzw_tf_t = Tzw_0;
        for t = 1:N_tf
            Tzw_t = sys.C1 * Rvar(:, t*n+1:(t+1)*n) * sys.B1 + ...
                    sys.D12 * Mvar(:, t*n+1:(t+1)*n) * sys.B1 + ...
                    sys.C1 * Nvar(:, t*p+1:(t+1)*p) * sys.D21 + ...
                    sys.D12 * Lvar(:, t*p+1:(t+1)*p) * sys.D21;
            Tzw_tf_t = Tzw_tf_t + z_i^(-t) * Tzw_t;
        end
        LMI_t = [lambda * eye(sys.n_z), Tzw_tf_t; Tzw_tf_t', lambda * eye(sys.n_w)];
        constraints = [constraints, LMI_t >= 0];
    end
elseif strcmp(problem, 'hinf_old_method')
    % OLD Method, using the max_svd approximation. Deprecated.
    Tzw_0 = sys.D12 * Lvar(:, 1:p) * sys.D21 + sys.D11;
    Tzw = Tzw_0;
    for t = 1:N_tf
        Tzw_i = sys.C1 * Rvar(:, t*n+1:(t+1)*n) * sys.B1 + ...
                sys.D12 * Mvar(:, t*n+1:(t+1)*n) * sys.B1 + ...
                sys.C1 * Nvar(:, t*p+1:(t+1)*p) * sys.D21 + ...
                sys.D12 * Lvar(:, t*p+1:(t+1)*p) * sys.D21;
        Tzw = blkdiag(Tzw, Tzw_i);
    end
    cost = N_tf * norm(Tzw, 2);
elseif strcmp(problem, 'l1')
    Tzw_0 = sys.D12 * Lvar(:,1:p) * sys.D21 + sys.D11;

    M_matrix = Tzw_0;
    
    for t= 1:N_tf
        Tzw_t = sys.C1  * Rvar(:, t*n+1:(t+1)*n) * sys.B1 + ...
                sys.D12 * Mvar(:, t*n+1:(t+1)*n) * sys.B1 + ...
                sys.C1  * Nvar(:, t*p+1:(t+1)*p) * sys.D21 + ...
                sys.D12 * Lvar(:, t*p+1:(t+1)*p) * sys.D21;
        M_matrix = [M_matrix, Tzw_t];
    end
    cost = norm(M_matrix,'inf');
    dummy_var = sdpvar(1);
    cost = cost +dummy_var^2; %this is forced to use IPM instead of Simplex (which is very slow in this case). Still figure out a better way.
    
elseif strcmp(problem, 'l1_regret')
    error("L1 Regret is not implemented yet.");
else
    error("Problem is not supported. Please choose 'l1', 'hinf', 'h2' or 'hinf_old_method'.");
end