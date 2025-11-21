function [K, optimization_parameter, objective] = calculus_distributed(sys, problem, delays_matrix, options)
clear yalmip; %removing old possible yalmip variables.
problem = lower(problem);  %formatting type variable
method = lower(options.method);

if strcmp(method, 'youla')
    if strcmp(problem, 'h2')
        [optimization_parameter, objective] = calculus_distributed_H2_youla(sys,delays_matrix,options);
    elseif strcmp(problem, 'hinf')
        [optimization_parameter, objective] = calculus_distributed_Hinf_youla(sys,delays_matrix,options);
    else
        error("\nInserted problem to solve wrongly!");
    end

    K = get_K_given_Q(sys,optimization_parameter);

elseif strcmp(method, 'sls')
    % if there is the flag distributed_optimization set to true in options, then use the distributed optimization
    if isfield(options,'distributed_optimization') && options.distributed_optimization
        assert(strcmp(problem,'l1'), "Error! Distributed optimization available only for L1 optimization method!")
        % disp("USING DISTRIBUTED OPTIMIZATION FOR SLS SYNTHESIS!!!!")
        [K, optimization_parameter, objective, ~] = distributed_optimization_L1_SLS(sys, delays_matrix, options); %Notice an optional output is available to check evolutoin of the problem during iterations
    else
        % disp("USING CENTRALIZED OPTIMIZATION FOR SLS SYNTHESIS!!!!")
        [K, optimization_parameter, objective] = calculus_distributed_sls(sys,problem,delays_matrix,options);
    end


elseif strcmp(method, 'sampled_youla')
    disp("USING ORIGINAL METHOD WITH SAMPLED MATRICES. NO TFS!!!!")
    if strcmp(problem, 'h2')
        [K, optimization_parameter, objective] = calculus_distributed_H2_youla_sampled(sys,delays_matrix,options);
    elseif strcmp(problem, 'hinf')
        [K, optimization_parameter, objective] = calculus_distributed_Hinf_youla_sampled(sys,delays_matrix,options);
    else
        error("\nInserted problem to solve wrongly!");
    end

elseif strcmp(method, 'tf_sampled')
    disp("USING METHOD WITH SAMPLING TFs!!!!")
    if strcmp(problem, 'h2')
        [K, optimization_parameter, objective] = calculus_distributed_H2_youla_sampled_WITH_TF(sys,delays_matrix,options);
    elseif strcmp(problem, 'hinf')
        [K, optimization_parameter, objective] = calculus_distributed_Hinf_youla_sampled_WITH_TF(sys,delays_matrix,options);
    else
        error("\nInserted problem to solve wrongly!");
    end
else
    error("\nInserted method wrongly!")
end
end