function violations = check_sls_rowwise(Phi, sys, N_tf, tol)
%CHECK_SLS_ROWWISE  Verify sls_achievability_constraints_only_rowwise on Phi
%
%   violations = check_sls_rowwise(Phi, sys, N_tf, tol)
%
%   Inputs:
%     Phi   : full decision matrix of size (n+m)x((n+p)*(N_tf+1))
%     sys   : system struct with fields n,m,p,A,C2
%     N_tf  : horizon used in the constructor
%     tol   : numeric tolerance for equality
%
%   Output:
%     violations : struct of flags (true if constraint is violated)

if nargin<4; tol = 1e-6; end

n = sys.n; m = sys.m; p = sys.p;
A = sys.A; C = sys.C2;

% Split Phi into blocks
R = Phi(1:n,           1 :   n*(N_tf+1));
M = Phi(n+1:end,       1 :   n*(N_tf+1));
N = Phi(1:n,           n*(N_tf+1)+1 :   end);
L = Phi(n+1:end,       n*(N_tf+1)+1 :   end);

violations = struct();

% R(0)=0 and M(0)=0
violations.R0 = any(abs(R(:,1:n)) > tol,'all');
violations.M0 = any(abs(M(:,1:n)) > tol,'all');

% R(1)=N(0)*C + I,  M(1)=L(0)*C
R1_expected = N(:,1:p)*C + eye(n);
M1_expected = L(:,1:p)*C;
violations.R1 = any(abs(R(:,n+1:2*n) - R1_expected) > tol,'all');
violations.M1 = any(abs(M(:,n+1:2*n) - M1_expected) > tol,'all');

% For t = 2…N_tf
violations.Rt = false(N_tf-1,1);
violations.Mt = false(N_tf-1,1);
for t = 2:N_tf
    idxR_t   = t*n+1 : (t+1)*n;
    idxR_t_1 = (t-1)*n+1 : t*n;
    idxN_t_1 = (t-1)*p+1 : t*p;
    Rt_expected = R(:,idxR_t_1)*A + N(:,idxN_t_1)*C;
    Mt_expected = M(:,idxR_t_1)*A + L(:,idxN_t_1)*C;
    violations.Rt(t-1) = any(abs(R(:,idxR_t)   - Rt_expected) > tol,'all');
    violations.Mt(t-1) = any(abs(M(:,idxR_t)   - Mt_expected) > tol,'all');
end

% Last block: R_end*A + N_end*C == 0,  M_end*A + L_end*C == 0
idxR_end = N_tf*n+1 : (N_tf+1)*n;
idxN_end = N_tf*p+1 : (N_tf+1)*p;
violations.Rend = any(abs(R(:,idxR_end)*A + N(:,idxN_end)*C) > tol,'all');
violations.Mend = any(abs(M(:,idxR_end)*A + L(:,idxN_end)*C) > tol,'all');

% Print summary
fprintf('Constraint violations:\n');
fprintf('  R(0)=0  : %d\n', violations.R0);
fprintf('  M(0)=0  : %d\n', violations.M0);
fprintf('  R(1)=…  : %d\n', violations.R1);
fprintf('  M(1)=…  : %d\n', violations.M1);
for t = 2:N_tf
  fprintf('  R(t=%d) : %d    M(t=%d) : %d\n', ...
          t, violations.Rt(t-1), t, violations.Mt(t-1));
end
fprintf('  Rend     : %d\n', violations.Rend);
fprintf('  Mend     : %d\n', violations.Mend);

end