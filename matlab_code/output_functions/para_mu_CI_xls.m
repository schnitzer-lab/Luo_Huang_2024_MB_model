function xls_mat = para_mu_CI_xls(para_mu,para_CI,mat_lu_cell)
%{
This function converts the parameters into the cell for xls output.
%}
para_mat_cell = [parameter_vec2mat(para_mu,mat_lu_cell), ...
                 parameter_vec2mat(para_CI(:,1),mat_lu_cell), ...
                 parameter_vec2mat(para_CI(:,2),mat_lu_cell)];

xls_mat = nan(23,3);
%{
Anti-Hebbian amplitude of KCs and MBONs
Anti-Hebbian amplitude on PPL1-α3

KC to PPL1-γ1pedc (Attractive odor)
KC to PPL1-α’2α2 (Attractive odor)
KC to PPL1-α3 (Attractive odor)
KC to PPL1-γ1pedc (Repulsive odor)
KC to PPL1-α’2α2 (Repulsive odor)
KC to PPL1-α3 (Repulsive odor)

Initial value of KC to MBON-γ1pedc>α/β
Initial value of KC to MBON-α’2α2
Initial value of KC to MBON-α3 

MBON-γ1pedc>α/β to PPL1-γ1pedc
MBON-γ1pedc>α/β to PPL1-α’2α2
MBON-γ1pedc>α/β to PPL1-α3
MBON-γ1pedc>α/β to MBON-α2sc
MBON-γ1pedc>α/β to MBON-α3 

MBON-α2sc to PPL1-α’2α2
MBON-α2sc to PPL1-α3
MBON-α3 to PPL1-α3

τMBON-STM (KC to MBON-γ1pedc>α/β) 
τMBON-STM (KC to MBON-α’2α2 and KC to MBON-α3)
τMBON-LTM (KC to MBON-α’2α2 and KC to MBON-α3)
τKC-recover
%}
for col_i = 1:size(xls_mat,2)
    xls_mat(:,col_i) = ...
        [para_mat_cell{2,col_i}/3;% The parameter in the fitting is the AAH of 3 training pulses.
         para_mat_cell{3,col_i}(1,1) * 11.38/3;% The parameter in the fitting is the AAH of 3 training pulses.
         para_mat_cell{1,col_i}([1:3,7:9,4:6])';
         para_mat_cell{4,col_i}(4,[1:3 5 6])';
         para_mat_cell{4,col_i}(5,[2 3])';
         para_mat_cell{4,col_i}(6,3);
         para_mat_cell{5,col_i}';
         para_mat_cell{6,col_i};
        ];
end

xls_mat(xls_mat==0) = nan;

for row_i = 1:size(xls_mat,1)
    if xls_mat(row_i,2) > xls_mat(row_i,3)
        xls_mat(row_i,[2 3]) = xls_mat(row_i,[3 2]);
    end
end
end