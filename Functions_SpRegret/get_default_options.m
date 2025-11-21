function options = get_default_options(varargin)
% GET_DEFAULT_OPTIONS Returns a structure with default options for spatial regret calculations.
%
%   OPTIONS = GET_DEFAULT_OPTIONS() returns a structure with default options.
%
%   OPTIONS = GET_DEFAULT_OPTIONS('PARAM1', VAL1, 'PARAM2', VAL2, ...) specifies
%   optional parameter name/value pairs to override the default values.
%
%   Parameters:
%       'method'        - (string) Method to be used. Default is 'sampled_youla'.
%       'number_points' - (integer) Number of sampling points. Default is 10.
%       'N_tf'          - (integer) N_tf value. Default is 10.
%       'trim_tol'      - (double) Trim tolerance. Default is 1e-9.
%       'verbose'       - (integer) Verbosity level. Default is 0.
%       'solver'        - (string) Solver to be used. Default is 'mosek'.
%       'N_FIR_closed_loop' - (integer) N_FIR_closed_loop value. Default is 10.
%       'warmstart'     - (boolean) Warmstart option. Default is false.
%       'warmStartSolution' - (struct) WarmStartSolution structure. Default is struct('R', [], 'M', [], 'N', [], 'L', []).
%       'distributed_optimization' - (boolean) Flag for distributed optimization. Default is false.

%   Returns:
%       options - A structure containing the specified options and default values
%                 for optimization options.
%
%   Example:
%       opts = get_default_options('method', 'new_method', 'number_points', 20);
%       This will return a structure with 'method' set to 'new_method' and
%       'number_points' set to 20, while other options will have their default values.


    % Default values
    default_method = 'sampled_youla';
    default_number_points = 100;  % Default value from NUMBER_OF_SAMPLING_POINTS
    default_N_tf = 10;           % Default value from N_tf
    default_trim_tol = 1e-10;
    default_verbose = 0;
    default_solver = 'mosek';
    default_N_FIR_closed_loop = 10;  % default value
    default_warmstart = false;  % default value for warmstart
    default_warmStartSolution = struct('R', [], 'M', [], 'N', [], 'L', []);  % new default value for warmStartSolution
    default_distributed_optimization = false;  % default value for distributed_optimization

    % Create input parser
    p = inputParser;
    addParameter(p, 'method', default_method);
    addParameter(p, 'number_points', default_number_points);
    addParameter(p, 'N_tf', default_N_tf);
    addParameter(p, 'trim_tol', default_trim_tol);
    addParameter(p, 'verbose', default_verbose);
    addParameter(p, 'solver', default_solver);
    addParameter(p, 'N_FIR_closed_loop', default_N_FIR_closed_loop);  
    addParameter(p, 'warmstart', default_warmstart);  
    addParameter(p, 'warmStartSolution', default_warmStartSolution);  
    addParameter(p, 'distributed_optimization', default_distributed_optimization); 

    % Parse inputs
    parse(p, varargin{:});

    % Create options object
    options = sdpsettings('verbose', p.Results.verbose, 'solver', p.Results.solver);
    options.method = p.Results.method;
    options.number_points = p.Results.number_points;
    options.N_tf = p.Results.N_tf;
    options.trim_tol = p.Results.trim_tol;
    options.N_FIR_closed_loop = p.Results.N_FIR_closed_loop;  
    options.warmstart = p.Results.warmstart;  
    options.warmStartSolution = p.Results.warmStartSolution;  
    options.distributed_optimization = p.Results.distributed_optimization;  
    % options.mosek.MSK_DPAR_INTPNT_TOL_REL_GAP = 1e-8;
    % options.mosek.MSK_DPAR_INTPNT_TOL_DFEAS = 1e-8;
    % options.mosek.MSK_DPAR_INTPNT_TOL_PFEAS = 1e-8;
    % options.mosek.MSK_IPAR_INTPNT_MAX_ITERATIONS = 1000;
    % options = sdpsettings(options, ...
    %                        'mosek.MSK_IPAR_LOG', 10, ...
    %                        'mosek.MSK_IPAR_INTPNT_BASIS', 'MSK_ON');
end