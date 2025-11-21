% spregnorm - Computes the spatial regret norm of the system_LFT with respect to the oracle_LFT.
%
% Syntax: 
%   lambda_opt = spregnorm(system_LFT, oracle_LFT, number_points)
%
% Inputs:
%   system_LFT - The system LFT (Linear Fractional Transformation) matrix. A tf object.
%   oracle_LFT - The oracle LFT matrix. A tf object.
%   number_points (optional) - The number of points to sample in the unit circle. Default is 100.
%
% Outputs:
%   lambda_opt - The optimal value of lambda, which represents the spatial regret norm.
%
% Example:
%   system_LFT = [1 2; 3 4];
%   oracle_LFT = [5 6; 7 8];
%   number_points = 50;
%   lambda_opt = spregnorm(system_LFT, oracle_LFT, number_points);
%
% See also: sdpvar, evalfr, optimize
function lambda_opt = spregnorm(system_LFT, oracle_LFT, number_points)
    % This function computes the spatial regret norm of the system_LFT with respect to the oracle_LFT.

    clear yalmip;
    if nargin < 3
        number_points = 100;
    end
    [n_z, n_w] = size(system_LFT);

    % Define the optimization problem
    lambda = sdpvar(1, 1);
    objective = lambda;

    % Define the constraints
    constraints = [lambda>=0];

    % Samples the points in the unit circle
    z_vector = exp(1i*2*pi*(0:number_points-1)/number_points);

    trim_tolerance = 1e-10;
    
    for i=1:number_points
        z_i = z_vector(i);
    
        Oracle_tf_t = evalfr(oracle_LFT,z_i);
        system_LFT_t = evalfr(system_LFT,z_i);
        Delta_tf_t = Oracle_tf_t'*Oracle_tf_t;
        
        %Trimming matrices
        system_LFT_t = trim_matrix(system_LFT_t,trim_tolerance);
        Delta_tf_t = trim_matrix(Delta_tf_t,trim_tolerance);
        
        % LMI_t = lambda*eye(n_w) - system_LFT_t'*system_LFT_t + Oracle_tf_t'*Oracle_tf_t;
        LMI_t = [1*eye(n_z), system_LFT_t; system_LFT_t', lambda*eye(n_w) + Delta_tf_t];
        constraints = [constraints, LMI_t >= 0];
    end

    % Solve the optimization problem
    options = sdpsettings('verbose', 0, 'solver', 'mosek');
    sol = optimize(constraints, objective, options);
    if ~(sol.problem == 0)
        error("Error solving the optimization problem. Problem status: " +num2str(sol.problem) + " ( " + yalmiperror(sol.problem) + ")");
    end
    lambda_opt =value(lambda);
end