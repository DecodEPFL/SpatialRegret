function [A_trim] = trim_matrix(A, tol)
% TRIM_MATRIX - Removes small entries from a matrix, setting them to zero.
%
%   A_trim = trim_matrix(A, tol)
%
%   This function removes small entries (values below a specified tolerance) from a matrix A,
%   setting them to zero. The resulting trimmed matrix is returned as A_trim.
%
%   Input arguments:
%       - A: Input matrix to be trimmed.
%       - tol: (Optional) Tolerance value for considering entries as small. Default is 1e-5.
%
%   Output:
%       - A_trim: Matrix A after removing small entries below the specified tolerance.
%
%   Example:
%       A = [0.001, 0.02; 0.03, 0.0001];
%       A_trim = trim_matrix(A, 0.01);
%       % A_trim = [0, 0.02; 0.03, 0]
%
%   See also ABS.

if nargin<2
    tol = 1e-5;
end
    indices_zeros = find(abs(A)<tol);
    A(indices_zeros) = 0;
    A_trim = A;
end

