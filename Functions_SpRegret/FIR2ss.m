function system = FIR2ss(FIR_matrix,Ts)
% FIR2ss - Given a mxnxN matrix of the coefficients of a FIR tf, it generates a state-space representation of the system.
%
%   The FIR matrix represents the coefficients of a FIR filter for multiple inputs
%   and outputs. The resulting state-space system is returned as 'system'.
%
%   Inputs:
%       - FIR_matrix: A 3D matrix of size (p x m x N_plus_1), where p is the number
%         of outputs, m is the number of inputs, and N_plus_1 is the number of
%         filter coefficients plus 1.
%       - Ts: The sampling time of the system. If not provided, the default value
%         is -1.
%
%   Output:
%       - system: The resulting state-space system.
%
%   Example:
%       FIR_matrix = rand(2, 3, 5); % Random FIR matrix
%       Ts = 0.1; % Sampling time
%       system = FIR2ss(FIR_matrix, Ts); % Convert FIR matrix to state-space system
%
%   See also ss
%
if nargin<2
    Ts = -1;
end
[p,m,N_plus_1] = size(FIR_matrix);

N = N_plus_1 - 1;
if (N == 0)
    %The system is static so A,B,C == 0
    Dtot = FIR_matrix(:,:,1);
    Atot = 0;
    Btot = zeros(1,m);
    Ctot = zeros(p,1);
else
    %The system is not static
    Aij = diag(ones(N-1,1),1); % The A matrix for each input j output i!
    Atot = kron(eye(m), Aij);

    Bij = zeros(N,1); Bij(end,1) = 1; %Bij = [0 ... 0 1]'
    Btot = kron(eye(m), Bij);

    Ctot = zeros(p,m*N);
    Dtot = zeros(p,m);

    for i=1:p
        for j = 1:m
            col_start = (j-1)*N + 1;
            col_end = j*N;
            Ctot(i,col_start:col_end) = fliplr(reshape(FIR_matrix(i,j,2:end), [1, N]));
            Dtot(i,j) = FIR_matrix(i,j,1);
        end
    end
end
system = ss(Atot,Btot,Ctot,Dtot,Ts);
end