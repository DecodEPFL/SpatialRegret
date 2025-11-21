function G = random_mesh_grid(shape, p)
% Generates an m-by-n random mesh grid where each point may be connected to
% its closest neighbors with probability p.
%
% Inputs:
%   n - size of the grid (n x n)
%   p - probability of connection to each neighbor [0, 1]
%
% Output:
%   G - adjacency matrix representing the grid

% Number of nodes in the grid
assert(length(shape) == 2, "The input 'shape' should be a 2-D array, instead it has dimension= " + num2str(length(shape)));
m = shape(1);
n = shape(2);
num_nodes = m * n;

% Initialize adjacency matrix
G = zeros(num_nodes, num_nodes);

% Iterate over each node in the grid
for i = 1:m
    for j = 1:n
        % Current node index
        current_node = (i - 1) * n + j;

        % Neighboring nodes
        neighbors = [];
        if i > 1
            neighbors(end+1) = current_node - n; % Up
        end
        if i < m
            neighbors(end+1) = current_node + n; % Down
        end
        if j > 1
            neighbors(end+1) = current_node - 1; % Left
        end
        if j < n
            neighbors(end+1) = current_node + 1; % Right
        end

        % Randomly connect to neighbors based on probability p
        for neighbor = neighbors
            if rand <= p
                G(current_node, neighbor) = 1;
                G(neighbor, current_node) = 1; % Since the graph is undirected
            end
        end
    end
end
% %Promoting that each node has at least one connection
for i = 1:m
    for j = 1:n
        % Current node index
        current_node = (i - 1) * n + j;

        % Neighboring nodes
        neighbors = [];
        if i > 1
            neighbors(end+1) = current_node - n; % Up
        end
        if i < m
            neighbors(end+1) = current_node + n; % Down
        end
        if j > 1
            neighbors(end+1) = current_node - 1; % Left
        end
        if j < n
            neighbors(end+1) = current_node + 1; % Right
        end

        edges_node_ij = G(current_node, neighbors);

        if ~any(edges_node_ij)
            % Randomly connect one of the neighbors
            neighbor = neighbors(randi(length(neighbors)));
            G(current_node, neighbor) = 1;
            G(neighbor, current_node) = 1; % Since the graph is undirected
        end
    end
end
% %Guaranteeing that you have self-loops
G = bin(G + eye(num_nodes));


end