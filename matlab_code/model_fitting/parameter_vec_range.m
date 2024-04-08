function para_lu = parameter_vec_range(mat_lu_cell)
para_lu = zeros(0,2);
for mat_i = 1:numel(mat_lu_cell)
    is_diff_lu = mat_lu_cell{mat_i}(:,:,1) ~= mat_lu_cell{mat_i}(:,:,2); 
    para_lu_curr = zeros(sum(sum(is_diff_lu)),2);
    for col_i=1:2
        mat_curr = mat_lu_cell{mat_i}(:,:,col_i);
        para_lu_curr(:,col_i) = mat_curr(is_diff_lu);
    end
para_lu = [para_lu; para_lu_curr];    
end

is_reverse = para_lu(:,1)>para_lu(:,2);
para_lu(is_reverse,:) = para_lu(is_reverse,[2 1]);
end