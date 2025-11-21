function [K, Phis] = sls_postprocessing(R,M,N,L,N_tf,Ts)

z = tf('z',Ts);
[~,p1] = size(L);
[n,~] = size(R);
p      = p1/(N_tf+1);
m = size(M,1);

Rtf = zeros(n,n); Mtf = zeros(m,n); Ntf = zeros(n,p); Ltf = zeros(m,p);
% disp("    TODO: Change this with a more optimal algorithm one day.")
for k = 1:N_tf+1
    Rtf = Rtf + R(:,n*(k-1)+1:n*k)*z^(1-k);  
    Mtf = Mtf + M(:,n*(k-1)+1:n*k)*z^(1-k); 
    Ntf = Ntf + N(:,p*(k-1)+1:p*k)*z^(1-k);  
    Ltf = Ltf + L(:,p*(k-1)+1:p*k)*z^(1-k);
end

Phis.R = Rtf;    Phis.M = Mtf;    Phis.N = Ntf;    Phis.L = Ltf;



%% state-space realization for L, M, zR, zN
% L = ss(Zp,hatIp,hatL,L0)
% N = ss(Zp,hatIp,hatN,N0)
hatIp = zeros(N_tf*p,p);  hatIp(1:p,1:p) = eye(p);
Zp   = diag(ones(N_tf-1,1), -1);     % downshift operator
Zp   = kron(Zp,eye(p));
hatL = L(:,p+1:end);
L0   = L(:,1:p);
hatN = N(:,p+1:end);
N0   = N(:,1:p);                   % this one must be zero

% zM = ss(Zn,hatIn,hatM,M1)
% zR = ss(Zn,hatIn,hatR,R1)
hatIn  = zeros((N_tf-1)*n,n); hatIn(1:n,1:n) = eye(n);
Zn = diag(ones(N_tf-2,1), -1);     % downshift operator
Zn = kron(Zn,eye(n));
hatM   = M(:,2*n+1:end);
M1     = M(:,n+1:2*n);
R1     = R(:,n+1:2*n);           % this one must be identity
hatR   = R(:,2*n+1:end);


% state space realization
A = [Zn-hatIn*hatR -hatIn*hatN;
    zeros(p*N_tf,n*(N_tf-1))  Zp];
B = [zeros(n*(N_tf-1),p);hatIp];
C = [hatM-M1*hatR -M1*hatN+hatL];
D = L0;

K = ss(A,B,C,D,Ts);

end

