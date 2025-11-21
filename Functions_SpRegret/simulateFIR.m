function z = simulateFIR(G, w)
%SIMULATEFIR  Simulate an FIR MIMO system
%   z = SIMULATEFIR(G, w)
%   G: [n_z × n_w × N_FIR] impulse‐response tensor
%   w: [n_instants × n_w] input trajectory
%   z: [n_instants × n_z] output trajectory

[n_z,~,N_FIR] = size(G);
[n_instants,~] = size(w);

z = zeros(n_instants, n_z);
for t = 1:n_instants
    for k = 1:N_FIR
        idx = t - k + 1;
        if idx < 1, break; end
        % w(idx,:) is 1×n_w, G(:,:,k)' is n_w×n_z ⇒ 1×n_z
        z(t,:) = z(t,:) + w(idx,:) * squeeze(G(:,:,k))';
    end
end
end