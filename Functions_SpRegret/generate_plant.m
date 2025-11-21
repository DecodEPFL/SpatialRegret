function plant = generate_plant(Adjacency_matrix, T_s, is_state_feedback)

n_agents = size(Adjacency_matrix,1);

% Range of the random parameters of the system
kij_range= [10, 20];
mi_range = [0.1, 1];

percentage_di = 40/100;
min_crit_dumping  = 2*sqrt(kij_range(1)*mi_range(1));
max_crit_dumping  = 2*sqrt(kij_range(2)*mi_range(2));

di_range =  percentage_di* [min_crit_dumping, max_crit_dumping];
fprintf("Range of dampers: [%.4f , %.4f ] \n",di_range);

% Initialize parameters
kij = kij_range(1) + (kij_range(2) - kij_range(1)) * rand(n_agents);
kij = kij.*Adjacency_matrix;
di = di_range(1) + (di_range(2) - di_range(1)) * rand(n_agents, 1);
mi = mi_range(1) + (mi_range(2) - mi_range(1)) * rand(n_agents, 1);


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
    A_tot(2*i, 2*i-1) = -T_s * (sum(kij(i, :)) / mi(i));
    A_tot(2*i, 2*i) = 1 - T_s * (di(i) / mi(i));
    
    % B_tot(2*i, i) = T_s / mi(i);

    for j = 1:n_agents
        if Adjacency_matrix(i,j) && i~=j
            A_tot(2*i, 2*j-1) = T_s * (kij(i, j) / mi(i));
        end
    end
end
plant=ss(A_tot,B_tot,C_tot,D_tot,T_s);
end

