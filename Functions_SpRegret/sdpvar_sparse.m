function sparse_variable = sdpvar_sparse(mask_matrix)
    % Find the indices of the non-zero elements in mask
    [row, col] = find(mask_matrix);
    
    % Create an sdpvar vector with the number of variables equal to the number of non-zero elements in mask
    num_vars = nnz(mask_matrix);
    variables = sdpvar(num_vars, 1);
    
    % Create the sparse sdpvar matrix Cq with the sparsity pattern defined by mask
    sparse_variable = sparse(row, col, variables, size(mask_matrix, 1), size(mask_matrix, 2));

%     % % Second (more intuitive way): Generate a full matrix and then mask it using hadamard product
%     sparse_variable = sdpvar(size(mask_matrix,1), size(mask_matrix,2), 'full');
%     sparse_variable = sparse_variable.*mask_matrix;
end