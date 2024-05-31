clc,clear
file_parameter_list = ...
    {'Dx_steady_state_nonlinear_3_27-Mar-2023_2modules.mat'  
     'Dx_steady_state_nonlinear_3_27-Mar-2023_3modules.mat'};
parameter_folder = fullfile('..','..','data_and_parameters');
load(fullfile(parameter_folder,file_parameter_list{1}))
mat_lu_cell_2 = mat_lu_cell;
para_CI_2 = para_CI;
para_mu_2 = para_mu;
para_rand_2 = para_rand;

load(fullfile(parameter_folder,file_parameter_list{2}))
mat_lu_cell_3 = mat_lu_cell;
para_CI_3 = para_CI;
para_mu_3 = para_mu;
para_rand_3 = para_rand;

%%
para_ratio = nan(size(para_rand_3));
for rand_i = 1:size(para_rand_3,2)
    xls_cell_2_curr = para_mu_CI_xls(para_rand_2(:,rand_i),para_CI_2,mat_lu_cell_2);
    xls_cell_3_curr = para_mu_CI_xls(para_rand_3(:,rand_i),para_CI_3,mat_lu_cell_3);

    para_ratio(:,rand_i) = xls_cell_2_curr(:,1)./xls_cell_3_curr(:,1);
end
%%
%para_ratio_prctile = prctile(para_ratio,[25 50 75],2);
%para_ratio_prctile = prctile(para_ratio,[5 50 95],2);
para_ratio_prctile = prctile(para_ratio,[16 50 84],2);
para_ratio_prctile = para_ratio_prctile(~isnan(para_ratio_prctile(:,1)),:);
%{
para_name = {'A_{AH}(0)';
             'A_{AH}(\Deltat_{punish})w_{punish,3}';
             'w_{KD,1,1},w_{KD,2,1}';
             'w_{KD,1,3},w_{KD,2,3}';
             'w_{KD,3,1},w_{KD,4,1}';
             'w_{KD,3,3},w_{KD,4,3}';
             'w_{KM,initial,1}';
             'w_{KM,initial,3}';
             'w_{MD,1,1}';
             'w_{MD,1,3}';
             'w_{MM,1,3}';
             'w_{MD,3,3}';
             '\tau_{u,1}';
             '\tau_{u,2},\tau_{u,2}(t\leq3 hrs)';
             '\tau_{u,2},\tau_{u,2}(t>3 hrs)';
             '\tau_{KC,recover}'};
%}
para_name = {'{\itA_{AH}}(0)';
             '{\itA_{AH}}(\Delta{\itt_{punish}}){\itw}_{{\itpunish},3}';
             '{\itw}_{{\itKD},1,1},{\itw}_{{\itKD},2,1}';
             '{\itw}_{{\itKD},1,3},{\itw}_{{\itKD},2,3}';
             '{\itw}_{{\itKD},3,1},{\itw}_{{\itKD},4,1}';
             '{\itw}_{{\itKD},3,3},{\itw}_{{\itKD},4,3}';
             '{\itw}_{{\itKM},{\itinitial},1}';
             '{\itw}_{{\itKM},{\itinitial},3}';
             '{\itw}_{{\itMD},1,1}';
             '{\itw}_{{\itMD},1,3}';
             '{\itw}_{{\itMM},1,3}';
             '{\itw}_{{\itMD},3,3}';
             '\tau_{{\itu},1}';
             '\tau_{{\itu},2},\tau_{{\itu},2}({\itt} \leq 3 hrs)';
             '\tau_{{\itu},2},\tau_{{\itu},2}({\itt} > 3 hrs)';
             '\tau_{{\itKC},{\itrecover}}'};

x = (1:16)';
%% vertical version
%{
figure, hold on
%boxplot(para_ratio')

errorbar(x, para_ratio_prctile(:,2), ...
            para_ratio_prctile(:,2) - para_ratio_prctile(:,1),  ...
            para_ratio_prctile(:,3) - para_ratio_prctile(:,2), ...
            'x')
%plot(x, para_ratio_prctile)
plot([0 17],[1 1],'k--')

set(gca,'YScale','log','XTick',x,'XTickLabel',para_name,'XTickLabelRotation',90, ...
    'XLim',[0 17], 'YLim', [0.1 10])
ylabel('parameter ratio (2-module model/3-module model)')
%}
%% horizontal version
figure, hold on
%boxplot(para_ratio')
x = (1:16)';
errorbar(para_ratio_prctile(:,2), x, ...
            para_ratio_prctile(:,2) - para_ratio_prctile(:,1),  ...
            para_ratio_prctile(:,3) - para_ratio_prctile(:,2), ...
            'x','horizontal')
%plot(x, para_ratio_prctile)
plot([1 1],[0 17],'k--')

set(gca,'XScale','log','YTick',x,'YTickLabel',para_name, ...
    'XLim',[0.1 10],'YLim',[0 17],'YDir','reverse', ...
    'FontName','Times New Roman')
xlabel('parameter ratio (two-module model/three-module model)','fontname','Arial')

