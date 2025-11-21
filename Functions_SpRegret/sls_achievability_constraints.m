function [Const] = sls_achievability_constraints(vars, sys, N_tf)
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
    Const = [Const, (Mvar(:,1:n) == zeros(m,n)):'M(0) = 0'];
    Const = [Const, (Nvar(:,1:p) == zeros(n,p)):'N(0) = 0'];

    % second spectral component
    Const = [Const, (Rvar(:,n+1:2*n) == eye(n)):'R(1) = I'];
    Const = [Const, (Nvar(:,p+1:2*p) == B*Lvar(:,1:p)):'N(1) = B*L(0)'];
    Const = [Const, (Mvar(:,n+1:2*n) == Lvar(:,1:p)*C):'M(1) = L(0)*C'];

    % the rest of them
    for t = 2:N_tf
        Const = [Const, (Rvar(:,t*n+1:(t+1)*n) == A*Rvar(:,(t-1)*n+1:t*n) + B*Mvar(:,(t-1)*n+1:t*n)):['R(' num2str(t) ') = A R(' num2str(t-1) ') + B M(' num2str(t-1) ')']];
        Const = [Const, (Nvar(:,t*p+1:(t+1)*p) == A*Nvar(:,(t-1)*p+1:t*p) + B*Lvar(:,(t-1)*p+1:t*p)):['N(' num2str(t) ') = A N(' num2str(t-1) ') + B L(' num2str(t-1) ')']];
        % Const = [Const, (Rvar(:,t*n+1:(t+1)*n) == Rvar(:,(t-1)*n+1:t*n)*A + Nvar(:,(t-1)*p+1:t*p)*C):['R(' num2str(t+1) ') = R(' num2str(t) ')A + N(' num2str(t) ')C']];
        Const = [Const, (Mvar(:,t*n+1:(t+1)*n) == Mvar(:,(t-1)*n+1:t*n)*A + Lvar(:,(t-1)*p+1:t*p)*C):['M(' num2str(t) ') = M(' num2str(t-1) ')A + L(' num2str(t-1) ')C']];
    end

    % the last component
    % TODO: add possible tolerance here 
    Const = [Const, (A*Rvar(:,N_tf*n+1:end) + B*Mvar(:,N_tf*n+1:end) == 0):['A R(' num2str(N_tf) ') + B M(' num2str(N_tf) ') = 0']];
    Const = [Const, (A*Nvar(:,N_tf*p+1:end) + B*Lvar(:,N_tf*p+1:end) == 0):['[A N(' num2str(N_tf) ') + B L(' num2str(N_tf) ') = 0']];
    % Const = [Const, (Rvar(:,N_tf*n+1:end)*A + Nvar(:,N_tf*p+1:end)*C == 0):'R(end)A + N(end)C = 0'];
    Const = [Const, (Mvar(:,N_tf*n+1:end)*A + Lvar(:,N_tf*p+1:end)*C == 0):['M(' num2str(N_tf) ')A + L(' num2str(N_tf) ')C = 0']];
end









% function [] = sls_achievability_constraints(vars, sys, N_tf)
% % encoding achievability constraints in SLS using a single struct "vars"
%     % Extract variables from struct


%     n = sys.n;
%     m = sys.m;
%     p = sys.p;
%     A = sys.A;
%     B = sys.B2;
%     C = sys.C2;
%     % get Lvar (given) and determine the number of time‐steps
%     Lvar = vars.Lvar;
    
%     % Preallocate matrices defined as functions of Lvar
%     Rvar = zeros(n, n*(N_tf+1));
%     Mvar = zeros(m, n*(N_tf+1));
%     Nvar = zeros(n, p*(N_tf+1));
    
%     % Initial conditions (t = 0)
%     Rvar(:,1:n) = 0;
%     Mvar(:,1:n) = 0;
%     Nvar(:,1:p) = 0;
    
%     % t = 1
%     Rvar(:,n+1:2*n) = eye(n);
%     Mvar(:,n+1:2*n) = Lvar(:,1:p)*C;
%     Nvar(:,p+1:2*p) = B*Lvar(:,1:p);
    
%     % Recurrence for t = 2:N_tf
%     for t = 2:N_tf
%         R_prev = Rvar(:, (t-1)*n+1 : t*n);
%         M_prev = Mvar(:, (t-1)*n+1 : t*n);
%         N_prev = Nvar(:, (t-1)*p+1 : t*p);
%         L_prev = Lvar(:, (t-1)*p+1 : t*p);
        
%         Rvar(:, t*n+1:(t+1)*n) = R_prev*A + N_prev*C;
%         Mvar(:, t*n+1:(t+1)*n) = M_prev*A + L_prev*C;
%         Nvar(:, t*p+1:(t+1)*p) = A*N_prev + B*L_prev;
%     end

%     % Final conditions (t = 0)
% %     Const = [Const, ...
% %         (A*Rvar(:,N_tf*n+1:end) + B*Mvar(:,N_tf*n+1:end) == 0):'A R(end) + B M(end) = 0', ...
% %         (A*Nvar(:,N_tf*p+1:end) + B*Lvar(:,N_tf*p+1:end) == 0):'A N(end) + B L(end) = 0', ...
% %         (Rvar(:,N_tf*n+1:end)*A + Nvar(:,N_tf*p+1:end)*C == 0):'R(end)A + N(end)C = 0', ...
% %         (Mvar(:,N_tf*n+1:end)*A + Lvar(:,N_tf*p+1:end)*C == 0):'M(end)A + L(end)C = 0'];
    
%     % Final sanity check with tolerance tol
%     tol = 1e-8;
%     assert(norm(Rvar(:,1:n), 'fro') < tol, 'Sanity check failed: R(0) ~= 0');
%     assert(norm(Mvar(:,1:n), 'fro') < tol, 'Sanity check failed: M(0) ~= 0');
%     assert(norm(Nvar(:,1:p), 'fro') < tol, 'Sanity check failed: N(0) ~= 0');
%     assert(norm(Rvar(:,n+1:2*n) - eye(n), 'fro') < tol, 'Sanity check failed: R(1) ~= I');
%     assert(norm(Mvar(:,n+1:2*n) - Lvar(:,1:p)*C, 'fro') < tol, 'Sanity check failed: M(1) ~= L(0)*C');
%     assert(norm(Nvar(:,p+1:2*p) - B*Lvar(:,1:p), 'fro') < tol, 'Sanity check failed: N(1) ~= B*L(0)');
    
%     for t = 2:N_tf
%         R_prev = Rvar(:, (t-1)*n+1 : t*n);
%         M_prev = Mvar(:, (t-1)*n+1 : t*n);
%         N_prev = Nvar(:, (t-1)*p+1 : t*p);
%         L_prev = Lvar(:, (t-1)*p+1 : t*p);
%         R_cur  = Rvar(:, t*n+1:(t+1)*n);
%         M_cur  = Mvar(:, t*n+1:(t+1)*n);
%         N_cur  = Nvar(:, t*p+1:(t+1)*p);
        
%         chk = (norm(R_cur - (R_prev*A + N_prev*C), 'fro') < tol) && ...
%               (norm(R_cur - (A*R_prev + B*M_prev), 'fro') < tol) && ...
%               (norm(M_cur - (M_prev*A + L_prev*C), 'fro') < tol) && ...
%               (norm(N_cur - (A*N_prev + B*L_prev), 'fro') < tol);
%         assert(chk, ['Sanity check failed at time step t = ' num2str(t)]);
%     end
    
%     R_last = Rvar(:, end-n+1:end);
%     M_last = Mvar(:, end-n+1:end);
%     N_last = Nvar(:, end-p+1:end);
%     L_last = Lvar(:, end-p+1:end);
%     % Final constraint 1: Check A*R_last + B*M_last is close to 0
%     chk1 = (norm(A*R_last + B*M_last, 'fro') < tol);
%     assert(chk1, 'Sanity check failed on A*R_last + B*M_last');

%     % Final constraint 2: Check A*N_last + B*L_last is close to 0
%     chk2 = (norm(A*N_last + B*L_last, 'fro') < tol);
%     assert(chk2, 'Sanity check failed on A*N_last + B*L_last');

%     % Final constraint 3: Check R_last*A + N_last*C is close to 0
%     chk3 = (norm(R_last*A + N_last*C, 'fro') < tol);
%     assert(chk3, 'Sanity check failed on R_last*A + N_last*C');

%     % Final constraint 4: Check M_last*A + L_last*C is close to 0
%     chk4 = (norm(M_last*A + L_last*C, 'fro') < tol);
%     assert(chk4, 'Sanity check failed on M_last*A + L_last*C');

%     vars.Rvar = Rvar;
%     vars.Mvar = Mvar;
%     vars.Nvar = Nvar;

% end


