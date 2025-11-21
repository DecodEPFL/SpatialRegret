function lft_sys = my_lft(sys,K)

M_inv = (eye(sys.m) - K.D*sys.D22);
M_hat = eye(sys.p) + sys.D22/M_inv*K.D;

Atot = [sys.A + sys.B2/M_inv*K.D*sys.C2,       sys.B2/M_inv*K.C;
    K.B*M_hat*sys.C2,    K.A + K.B*sys.D22/M_inv*K.C];

Btot = [sys.B1+sys.B2/M_inv*K.D*sys.D21;
    K.B*M_hat*sys.D21];

Ctot = [sys.C1 + sys.D12/M_inv*K.D*sys.C2,     sys.D12/M_inv*K.C];

Dtot = sys.D11 + sys.D12/M_inv*K.D*sys.D21;

lft_sys = ss(Atot,Btot,Ctot,Dtot,sys.plant.Ts);
end

