function [F, L] = calculus_F_and_L(sys, adjacency_matrix, randomly)

if nargin ==2
    randomly = false;
end



%Checking system stability:
eigs_A  = eig(sys.A);
max_eig_A = max(abs(eigs_A));
if(max_eig_A < 1)
    disp("    System is already stable! (max_eig = " + num2str(max(abs(eigs_A))) + ")");
    F = zeros(sys.m,sys.n);
    L = zeros(sys.n,sys.p);

else
    fprintf("    System is unstable. Max eigenvalue of A: %.4f .\n",max_eig_A)
    fprintf("    Computing the F and L matrices.")
    if randomly
        N_MAX_ITERATIONS  = 1e7;
        F = get_F_randomly(sys,adjacency_matrix, N_MAX_ITERATIONS);
        fprintf(".");
        L = get_L_randomly(sys,adjacency_matrix, N_MAX_ITERATIONS);
        fprintf(".\n");
    else
        F = get_F_with_LMI(sys,adjacency_matrix);
        fprintf(".");
        L = get_L_with_LMI(sys,adjacency_matrix);
        fprintf(".\n");
    end

    eigs_with_F = abs(eig(sys.A + sys.B2*F));
    eigs_with_L = abs(eig(sys.A + L*sys.C2));
    assert(all(eigs_with_F < 1), "F is not stabilizing the system!");
    assert(all(eigs_with_L < 1), "L is not making the system detectable!");
    disp("    Final Eigs of the system:")
    fprintf("    >max_abs(eig(A+B2*F)): ")
    fprintf("%d \n",max(eigs_with_F));
    fprintf("    >max_abs(eig(A+L*C2)): ")
    fprintf("%d \n",max(eigs_with_L));


end

end



function [F] = get_F_randomly(sys,adjacency_matrix, N_MAX_ITERATIONS)

found_flag = false;
counter = 0;


% Computing the mask of the matrix F
N_agents = size(adjacency_matrix,1);
m_i = sys.m/N_agents;
p_i = sys.p/N_agents;
n_i = sys.n/N_agents;

mask = kron(adjacency_matrix,ones(m_i,n_i));


%% Stochastic Searching Routing
while(~found_flag && counter<N_MAX_ITERATIONS)
    F = randn(sys.m,sys.n).*mask;
    eigs_CL = eig(sys.A + sys.B2*F);
    if(max( abs(eigs_CL)) < 0.98)% && min(abs(eigs_CL))>=0.5 )
        %FOUND!
        found_flag = true; disp("Debug");
    end
    counter = counter + 1;

end
if ~found_flag
    error("A matrix F was not found in " + num2str(N_MAX_ITERATIONS) + " iterations :(")
end
end


function [L] = get_L_randomly(sys,adjacency_matrix, N_MAX_ITERATIONS)

found_flag = false;
counter = 0;


% Computing the mask of the matrix F
N_agents = size(adjacency_matrix,1);
m_i = sys.m/N_agents;
p_i = sys.p/N_agents;
n_i = sys.n/N_agents;

mask = kron(eye(N_agents),ones(n_i,p_i));


%% Stochastic Searching Routing
while(~found_flag && counter<N_MAX_ITERATIONS)
    L = randn(sys.n,sys.p).*mask;
    eigs_CL = eig(sys.A + L*sys.C2);
    if(max( abs(eigs_CL)) < 0.98)% && min(abs(eigs_CL))>=0.5 )
        %FOUND!
        found_flag = true; disp("Debug");
    end
    counter = counter + 1;

end
if ~found_flag
    error("A matrix L was not found in " + num2str(N_MAX_ITERATIONS) + " iterations :(")
end
end


function [F] = get_F_with_LMI(sys, adjacency_matrix)

% error("WARNING! COMPUTING F WITH LMI HAS NOT BEEN CODED YETF
% F dimension: m times n
n_agents = size(adjacency_matrix,1);
n = sys.n;
m = sys.m;
p = sys.p;

n_i = n/n_agents;
m_i = m/n_agents;
p_i = p/n_agents;

%% Optimization Variables
M = sdpvar(n,n, 'symmetric');
% G = sdpvar(sys.n, sys.n, 'full');
% R = sdpvar(m,n,'full');

G = [];
%% Constraints
constraints = [];
%%% G MUST BE BLOCK DIAGONAL
for i=1:n_agents
    G = blkdiag(G, sdpvar(n_i, n_i, 'full'));
end
%%% R MUST BE partitioned as A(G)
% Building the mask matrix for R
R_mask = zeros(m,n);
for i=1:n_agents
    for j=1:n_agents
        if adjacency_matrix(i,j)
            R_mask((i-1)*m_i+1:i*m_i, (j-1)*n_i+1:j*n_i) = ones(m_i, n_i);
        end
    end
end

R =sdpvar_sparse(R_mask);

% constraints = [constraints; R .* (~R_mask) == 0];
%LMIs
Matrix = [M, sys.A*G+ sys.B2*R;
    (sys.A*G+ sys.B2*R)', G+G'- M];

tol_LMI = 1e-4;
constraints = [constraints; Matrix>=tol_LMI*eye(size(Matrix,1))];

%% Objective
objective = 1;

%% Optimizing ...
options = sdpsettings('verbose', 0, 'solver', 'mosek');
sol = optimize(constraints, objective, options);
if ~(sol.problem == 0)
    error("Error during calculus of the prestabilizing matrix F. Problem status: " +num2str(sol.problem) + " ( " + yalmiperror(sol.problem) + ")");
end

R = value(R);
G = value(G);

F = R/G;
% F = trim_matrix(F);
end


function [L] = get_L_with_LMI(sys, adjacency_matrix)

% error("WARNING! COMPUTING L WITH LMI HAS NOT BEEN CODED YET!")
% L dimensions: n x p
n_agents = size(adjacency_matrix,1);
n = sys.n;
m = sys.m;
p = sys.p;

n_i = n/n_agents;
m_i = m/n_agents;
p_i = p/n_agents;

%% Optimization Variables
M = sdpvar(n, n, 'symmetric');
R = [];
G = [];         %%% G MUST BE BLOCK DIAGONAL
for i=1:n_agents
    G = blkdiag(G, sdpvar(n_i, n_i, 'symmetric'));
    R = blkdiag(R, sdpvar(p_i, n_i, 'full'));
end

%LMIs
Matrix = [M, sys.A'*G+ sys.C2'*R;
    (sys.A'*G+ sys.C2'*R)', G+G'-M];


constraints = Matrix >=1e-2*eye(size(Matrix,1));

%% Objective
% objective = norm(R,'fro');
objective = 1;

%% Optimizing ...
options = sdpsettings('verbose', 0, 'solver', 'mosek');
sol = optimize(constraints, objective, options);
if ~(sol.problem == 0)
    error("Error during calculus of the centralized. Problem status: " +num2str(sol.problem) + " ( " + yalmiperror(sol.problem) + ")");
end

R = value(R);
G = value(G);

L = (R/G)';
% L = trim_matrix(L);
end

