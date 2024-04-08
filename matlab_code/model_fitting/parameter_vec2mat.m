function para_mat_cell = parameter_vec2mat(para_vec,mat_lu_cell)
para_mat_cell = cell(size(mat_lu_cell));
index_end = 0;
for mat_i = 1:numel(mat_lu_cell)
    is_diff_lu = mat_lu_cell{mat_i}(:,:,1) ~= mat_lu_cell{mat_i}(:,:,2); 
    index_curr = index_end + (1:sum(sum(is_diff_lu)));
    index_end = index_curr(end);
    para_mat_cell{mat_i} = mat_lu_cell{mat_i}(:,:,1);
    para_mat_cell{mat_i}(is_diff_lu) = para_vec(index_curr);
end
end