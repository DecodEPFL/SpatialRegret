function plant = generate_plant_homogeneous(Adjacency_matrix, T_s, is_state_feedback)

n_agents = size(Adjacency_matrix,1);

% Range of the random parameters of the system
kij= 20;
di = 2;
mi = 1;

% Initialize parameters
kij = kij.*Adjacency_matrix;


% Initialize the overall system matrices
A_tot = zeros(2 * n_agents);
B_tot  = kron(eye(n_agents),[0;1]); %B_tot = zeros(2 * n_agents, n_agents);
if is_state_feedback
    C_tot = eye(2 * n_agents);
    D_tot = zeros(2 * n_agents, n_agents);
else
    C_tot = kron(eye(n_agents),[1,0]);
    D_tot = zeros(n_agents, n_agents);
end


% Fill A_tot
for i = 1:n_agents
    A_tot(2*i-1, 2*i-1) = 1;
    A_tot(2*i-1, 2*i) = T_s;
    A_tot(2*i, 2*i-1) = -T_s * (sum(kij(i, :)) / mi);
    A_tot(2*i, 2*i) = 1 - T_s * (di / mi);
    
    % B_tot(2*i, i) = T_s / mi(i);

    for j = 1:n_agents
        if Adjacency_matrix(i,j) && i~=j
            A_tot(2*i, 2*j-1) = T_s * (kij(i, j) / mi);
        end
    end
end
plant=ss(A_tot,B_tot,C_tot,D_tot,T_s);
end

