function A_z = FIR2tf(A,Ts, method)
% FIR2tf - Convert a matrix representation of a FIR filter to a transfer function object.
%
% Syntax:
%   A_z = FIR2tf(A, Ts, method)
%
% Inputs:
%   - A: Matrix representation of the FIR filter coefficients. The matrix
%        should have dimensions [n x n x N_tf_plus_1], where n is the number
%        of inputs/outputs and N_tf_plus_1 is the number of filter
%        coefficients plus 1. (TODO: case in which the filter has different number of inputs and outputs is not implemented yet)
%   - Ts: Sampling time of the filter. Default value is -1, which indicates
%         an unspecified sampling time.
%   - method: Method for converting the FIR filter to a transfer function.
%             Available options are:
%             - "fir" (default): if A is a FIR filter.
%             - "inverse": if A is the ctranspose of a FIR filter.
%             - "hermitian":  if A is a para-hermitian polynomial matrix.
%
% Output:
%   - A_z: Transfer function representation of the FIR filter.
%
% Notes:
%   - The function uses the MATLAB Control System Toolbox function 'tf' to
%     create the transfer function.
%   - The 'stackElements' function is used internally to stack the elements
%     of the matrix representation along the third dimension in a cell matrix.
%
% Example:
%   % Create a matrix representation of a FIR filter
%   A = zeros(2, 2, 5);
%   A(1, 1, 1) = 1;
%   A(2, 2, 3) = 2;
%   A(2, 1, 4) = 3;
%   A(1, 2, 5) = 4;
%
%   % Convert the FIR filter to a transfer function
%   A_z = FIR2tf(A);
%
%   % Display the transfer function
%   disp(A_z);
%
% See also: tf, stackElements
if nargin<3
    method = "fir";
end
if nargin<2
    Ts = -1;
end
method = lower(convertCharsToStrings(method));


switch method
    case "fir"
        cell_vers_A = stackElements(A);
        [m, n, N_tf_plus_1] = size(A);
        den = cell(n,n);
        for i = 1:m
            for j = 1:n
                num_i = cell_vers_A{i,j};
                index_to_simplify = 0;
                for k = 0 : length(num_i) - 1
                    if(num_i(end-k) == 0)
                        index_to_simplify = index_to_simplify + 1;
                    else
                        break;
                    end
                end
                cell_vers_A{i,j} = cell_vers_A{i,j}(1:N_tf_plus_1-index_to_simplify);

                den{i,j} = cell_vers_A{i,j}*0; den{i,j}(1) = 1;

            end
        end

        A_z = tf(cell_vers_A, den, Ts);


    case "inverse"
        %we need to flip the coefficients individually
        cell_vers_A = stackElements(A);
        [m, n, ~] = size(A);
        assert(m==n, "The matrix must be square. 'Inverse' method with non square matrices is not implemented yet.");
        for i=1:n
            for j=1:n
                cell_vers_A{i,j} = fliplr(cell_vers_A{i,j});
            end
        end
        den = cell(n,n); den(:) = {1};
        A_z = tf(cell_vers_A, den, Ts);


    case "hermitian"
        cell_vers_A = stackElements(A);
        [m, n , doubleN_plus_1] = size(A);
        N = (doubleN_plus_1 - 1)/2;
        den = cell(n,n);
        assert(m==n, "The matrix must be square. 'Hermitian' method with non square matrices is not implemented yet.");
        for i = 1:n
            for j= 1:n
                num_i = cell_vers_A{i,j};
                if num_i == 0 %means empty tf
                    cell_vers_A{i,j} = [0];
                    den{i,j} = [1];
                else
                    index_to_simplify = 0;
                    for k = 1:N
                        if (num_i(k) == 0)
                            index_to_simplify = index_to_simplify + 1;
                        else
                            break;
                        end
                    end
                    reversed_and_flipped_num_i = cell_vers_A{i,j}(1+index_to_simplify:end);
                    cell_vers_A{i,j} = fliplr(reversed_and_flipped_num_i);
                    grade_den = N+1-index_to_simplify;
                    den{i,j} = zeros(1,grade_den); den{i,j}(1,1)=1;
                end
            end
        end
        A_z = tf(cell_vers_A, den, Ts);

    otherwise
        error("Unknown Method. Inserted :" + method);
end


end




function cell_result = stackElements(A)
% stackElements - Stack elements of a 3D matrix along the third dimension and store the result in a cell array.
%
% Syntax:
%   cell_result = stackElements(A)
%
% Input Arguments:
%   - A: 3D matrix of size (n x n x m), where n is the size of the first two dimensions and m is the size of the third dimension.
%
% Output Arguments:
%   - cell_result: Cell array of size (n x n), where each cell contains either a scalar 0 or a column vector representing the stacked elements of A along the third dimension.
%
% Example:
%   A = rand(3, 3, 4); % Create a 3D matrix of size (3 x 3 x 4)
%   cell_result = stackElements(A); % Stack the elements of A along the third dimension and store the result in a cell array
%
%   % Display the result
%   for i = 1:size(cell_result, 1)
%       for j = 1:size(cell_result, 2)
%           disp(cell_result{i, j});
%       end
%   end

% Get the size of the matrix
[m, n, ~] = size(A);

% Initialize the result cell array
cell_result = cell(m, n);

% Iterate over each element in the first two dimensions
for i = 1:m
    for j = 1:n
        % Stack the elements along the third dimension
        stackedElements = squeeze(A(i, j, :));

        % Check if all elements are zero
        if all(stackedElements == 0)
            cell_result{i, j} = 0;
        else
            cell_result{i, j} = stackedElements';
        end
    end
end
end
