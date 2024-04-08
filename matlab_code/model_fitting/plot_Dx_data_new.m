function [error_bar_h, fit_h]=plot_Dx_data_new(axis_h,Dx_mean,Dx_SEM,Dx_DAN_MBON)%training_para
%{
para_mat_cell_all = parameter_vec2mat(para_vec,mat_lu_cell);
        para_mat_cell = para_mat_cell_all;
        para_mat_cell{1} = para_mat_cell{1}(:,1:6);
%{
if exist('training_para','var')
    if training_para.nonlinear
        
        Dx_DAN_MBON = Dx_steady_state_MBON_ReLu(para_mat_cell,training_para);
    else
        Dx_DAN_MBON = Dx_steady_state_linear_predict(para_mat_cell,training_para);
    end
else
    Dx_DAN_MBON = Dx_steady_state_linear(para_mat_cell);
end
        %}
        
Dx_DAN_MBON = Dx_fun(para_mat_cell);
%}
%%
title_list = ...
    {'DAN-\gamma1pedc', ...'PPL1-\gamma1pedc', ...
     'DAN-\alpha''2\alpha2', ...'PPL1-\alpha''2\alpha2', 
     'DAN-\alpha3'; ...'PPL1-\alpha3'; ...
     'MBON-\gamma1pedc>\alpha/\beta', ...
     'MBON-\alpha2sc', ...
     'MBON-\alpha3'};
error_bar_h = zeros([size(axis_h),2]);
fit_h = zeros([size(axis_h),2]);
line_color_list='rb';
for row_i=1:2
    for col_i=1:3
        error_bar_h(row_i,col_i,:) = errorbar(axis_h(row_i,col_i), ...
                 permute(Dx_mean((row_i-1)*3+col_i,:,:),[3 2 1]), ...
                 permute( Dx_SEM((row_i-1)*3+col_i,:,:),[3 2 1]));
        %errorbar(axis_h(row_i,col_i), ...
        %         Dx_mean((row_i-1)*3+col_i,2,:), ...
        %         Dx_SEM((row_i-1)*3+col_i,2,:),'b')
        fit_h(row_i,col_i,:) = plot(axis_h(row_i,col_i), ...
             permute(Dx_DAN_MBON((row_i-1)*3+col_i,:,:),[3 2 1]),'--x');
        for CS_i=1:2
            set(error_bar_h(row_i,col_i,CS_i),'Color',line_color_list(CS_i))
            set(fit_h(row_i,col_i,CS_i),'Color',line_color_list(CS_i))
        end
        title(axis_h(row_i,col_i),title_list(row_i,col_i))
        if col_i == 1
            ylabel(axis_h(row_i,col_i),'\Deltafiring rate')
        end
            xlabel(axis_h(row_i,col_i),'Imaging session #')
        %{
        %if row_i == 2
        set(axis_h(row_i,col_i),'XTick',1:4, ...
            'XTickLabel',{'Test 1','Test 2','Test 3','Test 4'},...{'pre-training','after 3x training','after 6x training','after 1 h'}, ...
            'XTickLabelRotation',45)
        %end
        %}
    end
end
set(axis_h(1,:),'YLim',[-8 8])
set(axis_h(2,:),'YLim',[-35 40])
set(axis_h,'XLim',[1 size(Dx_mean,3)])
end