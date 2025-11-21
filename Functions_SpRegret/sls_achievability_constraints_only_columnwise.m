function [Const] = sls_achievability_constraints_only_columnwise(vars, sys, N_tf)
% encoding achievability constraints in SLS using a single struct "vars"
    
    n = sys.n;
    m = sys.m;
    p = sys.p;
    A = sys.A;
    B = sys.B2;
    C = sys.C2;
    Const = [];
    % Extract variables from struct
    Lvar = vars.Lvar;
    Rvar = vars.Rvar;
    Nvar = vars.Nvar;
    Mvar = vars.Mvar;

    % first spectral component
    Const = [Const, (Rvar(:,1:n) == zeros(n)):'R(0) = 0'];
    % Const = [Const, (Mvar(:,1:n) == zeros(m,n)):'M(0) = 0'];
    Const = [Const, (Nvar(:,1:p) == zeros(n,p)):'N(0) = 0'];

    % second spectral component
    % Const = [Const, (Rvar(:,n+1:2*n) == eye(n)):'R(1) = I'];
    % Const = [Const, (Rvar(:,n+1:2*n) == Nvar(:,1:p)*sys.C2 + eye(n)):'R(1) = N(0)*C + I'];
    Const = [Const, (Rvar(:,n+1:2*n) == B*Mvar(:,1:n)+ eye(n)):'R(1) = B2*M(0) + I'];
    Const = [Const, (Nvar(:,p+1:2*p) == B*Lvar(:,1:p)):'N(1) = B*L(0)'];
    % Const = [Const, (Mvar(:,n+1:2*n) == Lvar(:,1:p)*C):'M(1) = L(0)*C'];

    % the rest of them
    for t = 2:N_tf
        Const = [Const, (Rvar(:,t*n+1:(t+1)*n) == A*Rvar(:,(t-1)*n+1:t*n) + B*Mvar(:,(t-1)*n+1:t*n)):['R(' num2str(t) ') = A R(' num2str(t-1) ') + B M(' num2str(t-1) ')']];
        Const = [Const, (Nvar(:,t*p+1:(t+1)*p) == A*Nvar(:,(t-1)*p+1:t*p) + B*Lvar(:,(t-1)*p+1:t*p)):['N(' num2str(t) ') = A N(' num2str(t-1) ') + B L(' num2str(t-1) ')']];
        % Const = [Const, (Rvar(:,t*n+1:(t+1)*n) == Rvar(:,(t-1)*n+1:t*n)*A + Nvar(:,(t-1)*p+1:t*p)*C):['R(' num2str(t+1) ') = R(' num2str(t) ')A + N(' num2str(t) ')C']];
        % Const = [Const, (Mvar(:,t*n+1:(t+1)*n) == Mvar(:,(t-1)*n+1:t*n)*A + Lvar(:,(t-1)*p+1:t*p)*C):['M(' num2str(t) ') = M(' num2str(t-1) ')A + L(' num2str(t-1) ')C']];
    end

    % the last component
    % TODO: add possible tolerance here 
    Const = [Const, (A*Rvar(:,N_tf*n+1:end) + B*Mvar(:,N_tf*n+1:end) == 0):['A R(' num2str(N_tf) ') + B M(' num2str(N_tf) ') = 0']];
    Const = [Const, (A*Nvar(:,N_tf*p+1:end) + B*Lvar(:,N_tf*p+1:end) == 0):['[A N(' num2str(N_tf) ') + B L(' num2str(N_tf) ') = 0']];
    % Const = [Const, (Rvar(:,N_tf*n+1:end)*A + Nvar(:,N_tf*p+1:end)*C == 0):'R(end)A + N(end)C = 0'];
    % Const = [Const, (Mvar(:,N_tf*n+1:end)*A + Lvar(:,N_tf*p+1:end)*C == 0):['M(' num2str(N_tf) ')A + L(' num2str(N_tf) ')C = 0']];
end