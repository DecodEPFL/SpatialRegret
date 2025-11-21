function plant = mass_spring_damper_system(n_agents,vals, Ts, is_attached)
if nargin<4
    is_attached = true;
end

A = zeros(n_agents*2, n_agents*2);

% First Car
if is_attached
    A(1:2, 1:2) = [0, 1; -2*(vals.K(1)/vals.M(1)), -(vals.Fr(1) + vals.C(1))/vals.M(1)]; % first block
else
    A(1:2, 1:2) = [0, 1; -(vals.K(1)/vals.M(1)), -(vals.Fr(1) + vals.C(1))/vals.M(1)]; % first block
end

A(1:2, 3:4) = [0, 0; (vals.K(1)/vals.M(1)), (vals.C(1)/vals.M(1))]; % next row block

for i = 2:n_agents
    if i ~= n_agents
        A(2*i-1:2*(i),2*(i-1)-1:2*(i-1)) = [0, 0; (vals.K(i-1)/vals.M(i)), (vals.C(i-1)/vals.M(i))]; % previous block
        A(2*i-1:2*(i),2*(i-1)+1:2*(i)) = [0, 1; -(vals.K(i-1) + vals.K(i))/vals.M(i), -(vals.Fr(i) + vals.C(i-1) + vals.C(i))/vals.M(i)]; % center block
        A(2*i-1:2*(i),2*(i-1)+3:2*(i-1)+4) = [0, 0; (vals.K(i)/vals.M(i)), (vals.C(i)/vals.M(i))]; % next block
    else
        % Leader/Last car
        A(2*i-1:2*(i),2*(i-1)-1:2*(i-1)) = [0, 0; (vals.K(i-1)/vals.M(i)), (vals.C(i-1)/vals.M(i))]; % previous row block
        A(2*i-1:2*(i),2*(i-1)+1:2*(i)) = [0, 1; -(vals.K(i-1)/vals.M(i)), -(vals.Fr(i) + vals.C(i-1))/vals.M(i)]; % final block
    end
end

B = zeros(2*n_agents, n_agents);
for i = 1:n_agents
    B(2*i-1:2*i, i) = [0; 1/vals.M(i)];
end
% B(1:2,1)=0;

is_state_feedback = false;

if is_state_feedback
    C = eye(2*n_agents);
    D = zeros(2*n_agents, n_agents);
else
    C = kron(eye(n_agents),[1,0]);
    D = zeros(n_agents, n_agents);
end
eigs_continouns = eig(A);

if max(real(eigs_continouns))>=0
    error("System is unstable even in CT");
end

A_dt = A*Ts + eye(size(A));
B_dt = B*Ts;
plant = ss(A_dt, B_dt, C, D,Ts);
end
