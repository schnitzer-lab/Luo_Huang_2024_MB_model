function [fig_h, error_bar_h, fit_h] = plot_Dx_data(Dx_mean,Dx_SEM,Dx_DAN_MBON_median_eb)
%{
This function plot the original data and fitting curve.
Dx_mean: the mean of original data
Dx_SEM: the SEM of original data
Dx_DAN_MBON_median_eb: the fitting curve. 
If the 5th dimension contains 3 elements, they are 50th, 16th, 84th
percentiles.
%}
fig_h = [0 0];
is_simplified_model = all(reshape(isnan(Dx_mean([2 5],:,:,:)),[],1));
for fig_i = 1:2
    fig_h(fig_i) = figure;
    axis_h = zeros(3,2);
    for axis_i = 1:numel(axis_h)
        axis_h(axis_i)=subplot(size(axis_h,2),size(axis_h,1),axis_i); 
        hold on;
    end
    axis_h=axis_h';
    [error_bar_h(:,:,:,fig_i), fit_h(:,:,:,fig_i)] = ...
        plot_Dx_data_new(axis_h, ...
        Dx_mean(:,:,:,fig_i),Dx_SEM(:,:,:,fig_i), ...
        Dx_DAN_MBON_median_eb(:,:,:,fig_i,:));
    if is_simplified_model
        delete(axis_h(:,2))
    end
end
end


function [error_bar_h, fit_h] = plot_Dx_data_new(axis_h,Dx_mean,Dx_SEM,Dx_DAN_MBON)%training_para
%{
This funciton plot the original data and the curve fitting data.
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

        if size(Dx_DAN_MBON,5) == 1
            fit_h(row_i,col_i,:) = plot(axis_h(row_i,col_i), ...
                permute(Dx_DAN_MBON((row_i-1)*3+col_i,:,:),[3 2 1]),'--x');
        else
            fit_h(row_i,col_i,:) = errorbar(axis_h(row_i,col_i), ...
                [1:size(Dx_DAN_MBON,3);1:size(Dx_DAN_MBON,3)]', ...
                permute(Dx_DAN_MBON((row_i-1)*3+col_i,:,:,:,1),[3 2 1]), ...
                permute(diff(Dx_DAN_MBON((row_i-1)*3+col_i,:,:,:,[2 1]),[],5),[3 2 1]), ...
                permute(diff(Dx_DAN_MBON((row_i-1)*3+col_i,:,:,:,[1 3]),[],5),[3 2 1]), ...
                '--x');
        end
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