function UB = lInfInducedNorm(G, N, L)
% lInfInducedNorm
%
%   UB = lInfInducedNorm(G, N, L) computes an
%   estimated upper bound on the l_inf induced norm of the
%   MIMO discrete-time LTI system G defined by:
%
%       x[k+1] = A*x[k] + B*u[k]
%       y[k]   = C*x[k] + D*u[k]
%
%   The l_inf induced norm is reinterpreted as the infinity norm
%   of an infinite Toeplitz matrix with block entries:
%       [ D         0      0     ... ]
%       [ C*B       D      0     ... ]
%       [ C*A*B     C*B    D     ... ]
%       [ C*A^2*B   C*A*B  C*B   ... ]
%       [  ...      ...    ...   ... ]
%
%   To approximate the infinite sum, we truncate the series at
%   a finite index N and then bound the tail (k = N+1 to infinity)
%   using an upper bound. In particular, for each output row i,
%
%       rho_minus(i) = sum(abs(D(i,:))) + sum_{k=0}^{N} sum(abs( (C*A^k*B)(i,:) ))
%
%   and an upper bound on the tail is given by
%
%       rho_plus(i)  = ( sum(abs( C(i,:) * A^(N+L) )) * norm(B, inf) ) / (1 - norm(A^L, inf))
%
%   where L is chosen such that norm(A^L, inf) < 1.
%
%   The estimated upper bound is then
%
%       UB = max_i { rho_minus(i) + rho_plus(i) }.
%
%   Inputs:
%       G          : system.
%       N          : truncation index for the finite part.
%       L          : tail parameter; should satisfy norm(A^L,inf) < 1.
%
%   Output:
%       UB         : estimated upper bound on ||G||_inf.
%
%   Example:
%       % Define a discrete-time system
%       A = [0.8 0.1; 0 0.9];
%       B = [1; 0.5];
%       C = [0.2 0.3];
%       D = 0.1;
%       G = ss(A,B,C,D,-1)
%       % Choose truncation parameter N and tail parameter L
%       N = 50;
%       L = 5;  % Make sure norm(A^L,inf) < 1
%       UB = lInfInducedNorm(G, N, L);
%

% Check dimensions
[A,B,C,D] = ssdata(G);
[n, m] = size(B);
[p, ~] = size(C);

% Preallocate finite part for each output row
rho_minus = zeros(p,1);

% Add the direct term D
for i = 1:p
    rho_minus(i) = rho_minus(i) + sum(abs(D(i,:)));
end

% Sum the contributions for k = 0 to N
for k = 0:N
    % Compute A^k
    if k == 0
        Ak = eye(n);
    else
        Ak = A^k;
    end
    % Compute the product C*A^k*B, size: p x m
    term = C * Ak * B;
    % For each output row i, add sum(abs(term(i,:)))
    for i = 1:p
        rho_minus(i) = rho_minus(i) + sum(abs(term(i,:)));
    end
end

% Now, compute the tail bound for each output row using parameter L.
% Check that norm(A^L, inf) < 1
A_L = A^L;
normA_L = norm(A_L, inf);
if normA_L >= 1
    error('Choose L such that norm(A^L, inf) < 1.');
end

% Compute norm(B,inf) (maximum absolute row sum)
normB = norm(B, inf);

rho_plus = zeros(p,1);
for i = 1:p
    % Compute the row vector: C(i,:) * A^(N+L)
    ANL = A^(N+L);
    tail_term = C(i,:) * ANL;
    % Its l1 norm is the sum of absolute values:
    norm_tail = sum(abs(tail_term));
    % Tail bound for row i:
    rho_plus(i) = (norm_tail * normB) / (1 - normA_L);
end

% The estimated upper bound is the maximum (over i) of (rho_minus(i) + rho_plus(i))
row_estimates = rho_minus + rho_plus;
UB = max(row_estimates);

end
