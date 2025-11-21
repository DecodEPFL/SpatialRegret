function K = get_K_given_Q(sys,Q)
Aj = sys.A + sys.B2*sys.F + sys.L*sys.C2 + sys.L*sys.D22*sys.F;
Bj = sys.B2+sys.L*sys.D22;
Cj = -sys.C2-sys.D22*sys.F;
Dj = -sys.D22;


M_K = inv(eye(sys.m) - Q.D*Dj);
M_K_hat = eye(sys.p) + Dj*M_K*Q.D;

Atot = [Aj + Bj*M_K*Q.D*Cj,       Bj*M_K*Q.C;
    Q.B*M_K_hat*Cj,    Q.A + Q.B*Dj*M_K*Q.C];

Btot = [-sys.L+Bj*M_K*Q.D;
    Q.B*M_K_hat];

Ctot = [sys.F+M_K*Q.D*Cj,  M_K*Q.C];

Dtot = M_K*Q.D;

Atot = trim_matrix(Atot); Btot = trim_matrix(Btot); Ctot = trim_matrix(Ctot); Dtot = trim_matrix(Dtot);

K = ss(Atot,Btot,Ctot,Dtot,sys.plant.Ts);
end

