function [Const] = sls_achievability_constraints_distrib(Phi_k, sys, N_tf_closed_loop, dim, dim_index_vector)
    n = sys.n;
    m = sys.m;
    p = sys.p;
    Const = [];
if strcmp(dim,"row")
    % Constraint : [R N][zI-A] = [I]
    %              [M L][-C2]    [0]
    F_matrix = [eye(n,n);
                zeros(m,n)];
    
    Phi_0_LEFT  = Phi_k(:, 1:n);
    Phi_0_RIGHT = Phi_k(:, (n*(N_tf_closed_loop+1)) + [1:p]);
    
    % Time t = 0
    Const = [Const, (Phi_0_LEFT == 0):' Phi_L(0) = 0'];

    % Time t = 1
    Phi_1_LEFT  = Phi_k(:, n+1:2*n);

    Const = [Const, (Phi_1_LEFT == Phi_0_RIGHT*sys.C2 + F_matrix(dim_index_vector(1):dim_index_vector(2),:)):'Phi_L(1) = Phi_R(0)*C2 + F(rj,:)'];
    % Time t = 2:N_tf_closed_loop

    for t = 2:N_tf_closed_loop
        Phi_t_LEFT = Phi_k(:, t*n+1:(t+1)*n);
        Phi_t_minus_1_LEFT = Phi_k(:, (t-1)*n+1:t*n);
        % Phi_t_RIGHT = Phi_k(:, (n*(N_tf_closed_loop+1)) + [t*p+1:(t+1)*p]);
        Phi_t_minus_1_RIGHT = Phi_k(:, (n*(N_tf_closed_loop+1)) + [(t-1)*p+1:t*p]);

        Const = [Const, (Phi_t_LEFT == Phi_t_minus_1_LEFT*sys.A + Phi_t_minus_1_RIGHT*sys.C2):'Phi_L(t) = Phi_L(t-1)*A + Phi_R(t-1)*C2'];
    end

    % Time t = N_tf_closed_loop
    Phi_N_tf_closed_loop_LEFT = Phi_k(:, N_tf_closed_loop*n+1:(N_tf_closed_loop+1)*n);
    Phi_N_tf_closed_loop_RIGHT = Phi_k(:, (n*(N_tf_closed_loop+1)) + [N_tf_closed_loop*p+1:(N_tf_closed_loop+1)*p]);
    Const = [Const, (Phi_N_tf_closed_loop_LEFT*sys.A == - Phi_N_tf_closed_loop_RIGHT*sys.C2):'Phi_L(N_tf_closed_loop)*A + Phi_R(N_tf_closed_loop)*C2 == 0'];


elseif strcmp(dim,"col") || strcmp(dim,"column")
    % Constraint : [zI-A  -B2][R N] = [I   0]
    %                         [M L]

    H_matrix = [eye(n) , zeros(n,p)];
    column_indices_left = dim_index_vector.left;
    column_indices_right = dim_index_vector.right;
    left_columns_indices = zeros(2, N_tf_closed_loop+1);
    right_columns_indices = zeros(2, N_tf_closed_loop+1);
    for t = 0:N_tf_closed_loop
        left_columns_indices(1,t+1) = column_indices_left(1) + t*n; 
        left_columns_indices(2,t+1) = column_indices_left(2) + t*n;
        right_columns_indices(1,t+1) = column_indices_right(1) + t*p;
        right_columns_indices(2,t+1) = column_indices_right(2) + t*p;
    end
    
    % Time t = 0
    column_indices_time0 = [left_columns_indices(1,1):left_columns_indices(2,1) , (n*(N_tf_closed_loop+1))+[right_columns_indices(1,1): right_columns_indices(2,1)]];
    Phi_0_UP = Phi_k(1:n, column_indices_time0);
    Phi_0_DOWN = Phi_k(n+1:end, column_indices_time0);
    
    Const = [Const, (Phi_0_UP == 0):'Phi_U(0) = 0'];
    
    % Time t = 1
    column_indices_time1 = [left_columns_indices(1,2):left_columns_indices(2,2) , (n*(N_tf_closed_loop+1))+[right_columns_indices(1,2): right_columns_indices(2,2)]];
    Phi_1_UP = Phi_k(1:n, column_indices_time1);
    % Phi_1_DOWN = Phi_k(n+1:end, column_indices_time1);
    H_matrix_col = H_matrix(:, dim_index_vector.indices_for_H_matrix(1):dim_index_vector.indices_for_H_matrix(2));

    Const = [Const, (Phi_1_UP == sys.B2*Phi_0_DOWN + H_matrix_col):'Phi_U(1) = B2*Phi_D(0) + H(:, cj)'];

    % Time t = 2:N_tf_closed_loop

    for t = 2:N_tf_closed_loop
        column_indices_time_t = [left_columns_indices(1,t+1):left_columns_indices(2,t+1) , (n*(N_tf_closed_loop+1))+[right_columns_indices(1,t+1): right_columns_indices(2,t+1)]];
        Phi_t_UP = Phi_k(1:n, column_indices_time_t);
        Phi_t_DOWN = Phi_k(n+1:end, column_indices_time_t);
        
        column_indices_time_t_minus_1 = [left_columns_indices(1,t):left_columns_indices(2,t) , (n*(N_tf_closed_loop+1))+[right_columns_indices(1,t): right_columns_indices(2,t)]];

        Phi_t_minus_1_UP = Phi_k(1:n, column_indices_time_t_minus_1);
        Phi_t_minus_1_DOWN = Phi_k(n+1:end, column_indices_time_t_minus_1);

        Const = [Const, (Phi_t_UP == sys.B2*Phi_t_minus_1_DOWN + sys.A*Phi_t_minus_1_UP):'Phi_U(t) = B2*Phi_D(t-1) + A*Phi_U(t-1)'];

    end
    % Time t = N_tf_closed_loop
    column_indices_time_N_tf_closed_loop = [left_columns_indices(1,N_tf_closed_loop+1):left_columns_indices(2,N_tf_closed_loop+1) , (n*(N_tf_closed_loop+1))+[right_columns_indices(1,N_tf_closed_loop+1): right_columns_indices(2,N_tf_closed_loop+1)]];
    Phi_N_tf_closed_loop_UP = Phi_k(1:n, column_indices_time_N_tf_closed_loop);
    Phi_N_tf_closed_loop_DOWN = Phi_k(n+1:end, column_indices_time_N_tf_closed_loop);

    Const = [Const, (sys.A*Phi_N_tf_closed_loop_UP + sys.B2*Phi_N_tf_closed_loop_DOWN == 0):'Phi_U(N_tf_closed_loop)*A + B2*Phi_D(N_tf_closed_loop) == 0'];

else
    error('Invalid dimension specified. Use "row" or "col".');
end