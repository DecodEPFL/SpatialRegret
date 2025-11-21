function plant = generate_tanks(Adjacency_matrix, T_s)

n_agents = size(Adjacency_matrix,1);

% Range of the random parameters of the system
ki_range = [0.80, 0.99];
Si_range = [0.5, 1.2];

% Initialize parameters
ki = ki_range(1) + (ki_range(2) - ki_range(1)) * rand(n_agents,1);
Si = Si_range(1) + (Si_range(2) - Si_range(1)) * rand(n_agents, 1);

Si_inv = 1./Si;
% Initialize the overall system matrices
A_tot = zeros(n_agents);
B_tot = diag(Si_inv);


% Fill A_tot
for i=1:n_agents
    A_tot(:,i) = ki(i,1);
end
% Changing the sign over the main diagonal
mask_sign = ones(n_agents,n_agents); mask_sign = mask_sign -2*eye(n_agents);
A_tot = A_tot.*mask_sign;

% Scaling everything by Si
Si_inv_matrix = kron(ones(1,n_agents),Si_inv);

A_tot = A_tot.*Si_inv_matrix;

% Removing the edges that do not exist
A_tot = A_tot.*Adjacency_matrix;

%matrix C_tot
C_tot = eye(n_agents);
D_tot = zeros(n_agents);



%%Discretizing using Forward Euler
Adt = A_tot*T_s + eye(n_agents);
Bdt = B_tot*T_s;
Cdt = C_tot;
Ddt = D_tot;

plant = ss(Adt,Bdt,Cdt,Ddt,T_s);
end