function final_L1norm = L1norm_FIR(sys, Phis, N_FIR, input_channel)
% L1norm_FIR computes the L1 norm (infinity norm) of the aggregate FIR transformation matrix.
%
% Syntax:
%   final_L1norm = L1norm_FIR(sys, Phis, N_FIR)
%   final_L1norm = L1norm_FIR(sys, Phis, N_FIR, input_channel)
%
% Description:
%   This function returns the infinity norm (L1 norm) of a FIR transfer function matrix G = [C1 D12]*[R N;M L]*[B1; D21] + D11.
%   The computation is based on the fact that the L1 norm of a system G with order N_FIR is equal to
%   ||G||_{L_1} = max_{i=1,...,p} \sum_{j=1,...,n} \sum_{k=0}^{N_FIR} |G_{ij}(k)|.
%
% Input Arguments:
%   sys           - Structure containing system matrices and dimensions:
%                     p: number of rows in related matrices.
%                     n: number of columns in related matrices.
%                     D12, D21, D11: system matrices.
%                     C1, B1: system matrices used in the transformation.
%
%   Phis          - Structure containing FIR matrices:
%                     L, R, M, N: matrices used to construct the transformation blocks.
%
%   N_FIR         - Scalar specifying the number of FIR terms (time steps) to consider.
%
%   input_channel - (Optional) Scalar indicating the specific input channel to select from
%                   Tzw_t. If input_channel > 0, only the corresponding column of Tzw_t is used.
%
% Output Arguments:
%   final_L1norm  - The infinity norm (L1 norm) of the final aggregated matrix M_matrix.
%
% Example:
%   % Define system structure with required matrices and dimensions
%   sys = struct('p', p_val, 'n', n_val, 'D12', D12, 'D21', D21, 'D11', D11, 'C1', C1, 'B1', B1);
%
%   % Define FIR matrices in Phis structure
%   Phis = struct('L', L_matrix, 'R', R_matrix, 'M', M_matrix, 'N', N_matrix);
%
%   % Set the number of FIR terms
%   N_FIR = 5;
%
%   % Compute the L1 norm of the FIR transformation matrix
%   L1_norm_value = L1norm_FIR(sys, Phis, N_FIR);

    p = sys.p;
    n = sys.n;
    % m = sys.m;
    if nargin < 4
        input_channel = 0;
    end

    Tzw_0 = sys.D12 * Phis.L(:,1:p) * sys.D21 + sys.D11;

    M_matrix = Tzw_0;
    
    for t= 1:N_FIR
        Tzw_t = sys.C1  * Phis.R(:, t*n+1:(t+1)*n) * sys.B1 + ...
                sys.D12 * Phis.M(:, t*n+1:(t+1)*n) * sys.B1 + ...
                sys.C1  * Phis.N(:, t*p+1:(t+1)*p) * sys.D21 + ...
                sys.D12 * Phis.L(:, t*p+1:(t+1)*p) * sys.D21;
        if input_channel > 0
            M_matrix = [M_matrix, Tzw_t(:, input_channel)];
        else
            M_matrix = [M_matrix, Tzw_t];
        end
    end
    final_L1norm = norm(M_matrix,'inf');
end
