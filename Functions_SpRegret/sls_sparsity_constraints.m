function L_sparsity = sls_sparsity_constraints(sys,delays_matrix,N_tf,N_tf_closed_loop)
% encoding achivability constraint in SLS
n = sys.n;    m = sys.m;    p = sys.p;
N_agents = size(delays_matrix,1);
p_i = p/N_agents;
m_i = m/N_agents;

% constraints = [];
L_sparsity = ones(m,p*(N_tf_closed_loop+1));

if (N_tf_closed_loop > N_tf)
    % Put zeros in the last N_tf_closed_loop - N_tf columns
    L_sparsity(: , p*(N_tf+1)+1 : end) = 0;
end

for i=1:N_agents
    for j=1:N_agents
        dist_ij = delays_matrix(i,j);

        if (dist_ij ~= 0 && dist_ij~= Inf && dist_ij<=N_tf) %you have a finite delay
            row_start = (i - 1) * m_i + 1; 
            row_end = i * m_i;
            for k = 1:dist_ij
                col_start = (k - 1) * p + (j - 1) * p_i + 1;
                col_end = (k - 1) * p + j * p_i;
                % fprintf("Finite Delay element(%d,%d) time %d --> (%d to %d,%d to %d) == 0 \n",i,j,k, row_start,row_end,col_start,col_end);
                % constraints = [constraints; Lvar(row_start:row_end,col_start:col_end) == 0];
                L_sparsity(row_start:row_end,col_start:col_end) = 0;
            end
        elseif (dist_ij == Inf || dist_ij > N_tf) %you have no communication!
            row_start = (i - 1) * m_i + 1;
            row_end = i * m_i;
            for k = 1:N_tf + 1
                col_start = (k - 1) * p + (j - 1) * p_i + 1;
                col_end = (k - 1) * p + j * p_i;
                % fprintf("Infinite Delay element(%d,%d) time %d --> (%d to %d,%d to %d) == 0 \n",i,j,k, row_start,row_end,col_start,col_end);
                % constraints = [constraints; Lvar(row_start:row_end,col_start:col_end) == 0];
                L_sparsity(row_start:row_end,col_start:col_end) = 0;
            end
        end

    end
end


end
