function Adjacency_matrix = generate_grid_16_adjacency()
% GENERATE_GRID_16_ADJACENCY - Creates the 16-agent grid adjacency matrix used in the paper
%
% This function generates the specific 4x4 grid topology used in Section IV.B
% of the paper. The structure includes custom edges for enhanced connectivity.
%
% Syntax:
%   Adjacency_matrix = generate_grid_16_adjacency()
%
% Outputs:
%   Adjacency_matrix - 16x16 binary adjacency matrix defining the network topology
%
% Example:
%   A = generate_grid_16_adjacency();
%   G = digraph(A ~= 0);
%   plot(G);
%
% See also: random_mesh_grid, generate_plant_homogeneous

% Define the adjacency matrix as specified in the paper
% This is a 16-node (4x4) grid with custom edge modifications
Adjacency_matrix = [
   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0
   0   1   1   0   0   1   0   0   0   0   0   0   0   0   0   0
   0   1   1   0   0   0   1   0   0   0   0   0   0   0   0   0
   0   0   0   1   0   0   0   1   0   0   0   0   0   0   0   0
   1   0   0   0   1   0   0   0   1   0   0   0   0   0   0   0
   0   1   0   0   0   1   0   0   0   1   0   0   0   0   0   0
   0   0   1   0   0   0   1   1   0   0   0   0   0   0   0   0
   0   0   0   1   0   0   1   1   0   0   0   0   0   0   0   0
   0   0   0   0   1   0   0   0   1   1   0   0   1   0   0   0
   0   0   0   0   0   1   0   0   1   1   1   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   1   1   1   0   0   1   0
   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0
   0   0   0   0   0   0   0   0   1   0   0   0   1   1   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   1   1   1   0
   0   0   0   0   0   0   0   0   0   0   1   0   0   1   1   1
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1
];

% Verify dimensions
assert(size(Adjacency_matrix, 1) == 16, 'Adjacency matrix must be 16x16');
assert(size(Adjacency_matrix, 2) == 16, 'Adjacency matrix must be 16x16');

end
