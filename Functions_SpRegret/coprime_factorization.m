function [P11,P12,P21] = coprime_factorization(sys, verbose)

if(nargin < 2)
    verbose  = 1;
end

Af= sys.A + sys.B2*sys.F;
Al= sys.A + sys.L*sys.C2;

% %% COMPUTING P11, P12, P21 --> [Aj,Bj,Cj,Dj]
% disp("NEW METHOD TO COMPUTE THE T transfer functions!!!!!!!")
% %P11
% A_tot = [sys.A+ sys.L*sys.C2, zeros(sys.n,sys.n);
%     -sys.L*sys.C2, sys.A + sys.B2*sys.F];
%
% P11.A = A_tot;
% P11.B = [sys.B1+sys.L*sys.D21; -sys.L*sys.D21];
% P11.C = [sys.C1 , sys.C1 + sys.D12*sys.F];
% P11.D = sys.D11;
% % P11.A = trim_matrix(P11.A); P11.B = trim_matrix(P11.B); P11.C = trim_matrix(P11.C); P11.D = trim_matrix(P11.D);
% P11 = ss(P11.A,P11.B,P11.C,P11.D,sys.plant.Ts);
% P11 = minreal(P11);
%
% %P12
% P12.A = A_tot;
% P12.B = [sys.B2; zeros(sys.n,sys.m)];
% P12.C = -[sys.C1 , sys.C1 + sys.D12*sys.F];
% P12.D = sys.D12;
% % P12.A = trim_matrix(P12.A); P12.B = trim_matrix(P12.B); P12.C = trim_matrix(P12.C); P12.D = trim_matrix(P12.D);
% P12 = ss(P12.A,P12.B,P12.C,P12.D,sys.plant.Ts);
% P12 = minreal(P12);
%
% %P21
% P21.A = A_tot;
% P21.B = [sys.B1+sys.L*sys.D21; -sys.L*sys.D21];
% P21.C = [sys.C2 , sys.C2 + sys.D22*sys.F];
% P21.D = sys.D21;
% % P21.A = trim_matrix(P21.A); P21.B = trim_matrix(P21.B); P21.C = trim_matrix(P21.C); P21.D = trim_matrix(P21.D);
% P21 = ss(P21.A,P21.B,P21.C,P21.D,sys.plant.Ts);
% P21 = minreal(P21);


%% COMPUTING P11, P12, P21 --> [Aj,Bj,Cj,Dj]
if( ~ isstable(sys.plant))
    %The system is not stable. A coprime factorization is required.
    %P11
    P11.A = [Af, -sys.B2*sys.F; zeros(size(Al,1), size(Af,2)), Al];
    P11.B = [sys.B1;sys.B1+sys.L*sys.D21];
    P11.C = [sys.C1 + sys.D12*sys.F, -sys.D12*sys.F];
    P11.D = sys.D11;
    % P11.A = trim_matrix(P11.A); P11.B = trim_matrix(P11.B); P11.C = trim_matrix(P11.C); P11.D = trim_matrix(P11.D);
    P11 = ss(P11.A,P11.B,P11.C,P11.D,sys.plant.Ts);

    %P12
    P12.A = Af;
    P12.B = sys.B2;
    P12.C = sys.C1 + sys.D12*sys.F;
    P12.D = sys.D12;
    % P12.A = trim_matrix(P12.A); P12.B = trim_matrix(P12.B); P12.C = trim_matrix(P12.C); P12.D = trim_matrix(P12.D);
    P12 = ss(P12.A,P12.B,P12.C,P12.D,sys.plant.Ts);

    %P21
    P21.A = Al;
    P21.B = sys.B1+sys.L*sys.D21;
    P21.C = sys.C2;
    P21.D = sys.D21;
    % P21.A = trim_matrix(P21.A); P21.B = trim_matrix(P21.B); P21.C = trim_matrix(P21.C); P21.D = trim_matrix(P21.D);
    P21 = ss(P21.A,P21.B,P21.C,P21.D,sys.plant.Ts);
else
    if verbose == 1
        disp("System is already stable. The matrices Ts are the plant itself!")
    end
    P11 = ss(sys.A,sys.B1,sys.C1,sys.D11,sys.plant.Ts);  % P11 ==  P11
    P12 = ss(sys.A,sys.B2,-sys.C1,-sys.D12,sys.plant.Ts); % P12 == -P12
    P21 = ss(sys.A,sys.B1, sys.C2,sys.D21,sys.plant.Ts); % P12 == P21
end
end