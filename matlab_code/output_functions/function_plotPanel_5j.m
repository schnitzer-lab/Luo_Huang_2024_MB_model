function function_plotPanel_5j(file_parameter,figure_name,figure_type)


%load('/Users/flyspike/Desktop/Fly leanring paper/fly model with CI/version2 - a2 a3 equal tau/figures/Dx_steady_state_nonlinear_3_18-Aug-2021')
%load('figures/Dx_steady_state_nonlinear_3_18-Aug-2021')
load(file_parameter)
addpath('customcolormap\')

%-----loading model fitted parameters
%para_mat_cell = parameter_vec2mat(x_lu,mat_lu_cell);
para_mat_cell_all = parameter_vec2mat(para_mu,mat_lu_cell);
switch figure_type
    case 'ACV_EtA'
        para_mat_cell = para_mat_cell_all;
        para_mat_cell{1} = para_mat_cell{1}(:,1:6);
    case 'OCT_BEN'
        para_mat_cell = para_mat_cell_all;
        para_mat_cell{1} = para_mat_cell{1}(:,[7:9 4:6]);
    otherwise
        error('figure_type should be ACV_EtA and OCT_BEN.')
end


%%
t_extinction_list = [nan 10 20 30 60 90 120]*60;
rest_t_list = (300:300:3600*3)-135;
output_data_matrix = nan(length(t_extinction_list),length(rest_t_list),3,6);
% the first two dimensions are row and column of the color maps, 
% the last two dimensions correspond to the row and column of subplots.

for t_extinction_i = 1:length(t_extinction_list)
    t_extinction_curr = t_extinction_list(t_extinction_i);
    if isnan(t_extinction_curr)
        % create the defualt experimental condition listexp_para
        [exp_para, exp_struct] = experimental_condition([]);
        
        exp_para.session_list = ...
            {'imaging';
            'training';
            'rest';
            'imaging'};
        
        exp_para.session(3).t_ISI = 135; % interval
        exp_para.session(3).n = 3; % training cycle 

        exp_para.rest_t = 0;
    else
        [exp_para, exp_struct] = experimental_condition([]);
        
        exp_para.session_list = ...
            {'imaging';
            'training';
            'rest';
            'extinction';
            'rest';
            'imaging'};
        
        exp_para.session(3).t_ISI = 135; % interval
        exp_para.session(3).n = 3; % training cycle 
        
        exp_para.session(4) = exp_para.session(3);
        exp_para.session(4).name = 'extinction';

        exp_para.rest_t = [t_extinction_curr; 0];
        
    end
    


    for rest_t_i = 1:length(rest_t_list)
        if isnan(t_extinction_curr)
            exp_para.rest_t = rest_t_list(rest_t_i);
        else
            exp_para.rest_t(end) = rest_t_list(rest_t_i) - exp_para.rest_t(1) ...
                - exp_para.session(4).n*(exp_para.session(4).t_CS_plus ...
                                        +exp_para.session(4).t_ISI ...
                                        +exp_para.session(4).t_CS_minus ...
                                        +exp_para.session(4).t_ISI);
        end
        
        [exp_para, exp_struct] = experimental_condition(exp_para);
        
        
        if exp_para.rest_t >= 0
            [Dx_DAN_MBON,data_all] = Dx_steady_state_MBON_0301_2023(para_mat_cell,exp_struct);
            % Dx_DAN_MBON: first dimension: DAN,MBON; second dimension: CS+, CS-
            output_data_matrix(t_extinction_i,rest_t_i,1:2,:) = permute(Dx_DAN_MBON(:,:,end),[4 3 2 1]);
        elseif rest_t_list(rest_t_i) <= t_extinction_curr
            output_data_matrix(t_extinction_i,rest_t_i,1:2,:) = ...
                output_data_matrix(1,rest_t_i,1:2,:);
        else
            output_data_matrix(t_extinction_i,rest_t_i,1:2,:) = nan;
        end
    end
    output_data_matrix(:,:,3,:) = output_data_matrix(:,:,1,:) - output_data_matrix(:,:,2,:);
    
end

%%
paper_size=[11 6];

fig_h=figure('Units','inches', ...
             'Papersize', paper_size, ...
             'PaperPosition', [0 0, paper_size], ...
             'Position',[0.5 0.5 paper_size]);

paper_margin=[(8.5-6.7)/2 0.5];
font_size=7;

% cMap = getPyPlot_cMap('RdBu_r');
J = customcolormap([0 0.25 0.375 0.5 0.625 0.75 1], {'#8B008B','#C71585','#FA8072','#FFFAFA','#48D1CC','#4682B4','#000000'});
% J = customcolormap([0 0.125 0.25 0.375 0.5 0.625 0.75 0.875 1], {'#8B008B','#C71585','#FA8072','#EEE8AA','#FFFAFA','#3CB371','#48D1CC','#4682B4','#000000'});

colormap(J);
axis_h = zeros(3,6);
for row_i = 1:3
    for col_i = 1:6
        axis_h(row_i,col_i) = subplot(3,6,(row_i-1)*6+col_i); hold on;
        if row_i <= 2
            imagesc(rest_t_list/3600,1:length(t_extinction_list), ...
                output_data_matrix(:,:,row_i,col_i)-Dx_DAN_MBON(6,row_i,1),[-40,40])
        else
            imagesc(rest_t_list/3600,1:length(t_extinction_list), ...
                output_data_matrix(:,:,row_i,col_i),[-40,40])
        end
    end
end


set(axis_h,'box','off','YDir','normal', ...
    'YLim',[0.5,length(t_extinction_list) + 0.5], ...
    'YTick',1:length(t_extinction_list), ...
    'YTickLabel',[{'Ctr'};cellstr(num2str(t_extinction_list(2:end)'/60))], ...
    'XTick',0:3)

cell_name = {'DAN-\gamma1', 'DAN-\alpha2', 'DAN-\alpha3', ...
             'MBON-\gamma1','MBON-\alpha2','MBON-\alpha3'};
for col_i = 1:size(axis_h,2)
    title_str = [cell_name{col_i},'(',figure_type(1:3),'/',figure_type(5:7),')'];
    title(axis_h(1,col_i),title_str)
end
%%
%{
if ~isempty(figure_name)
    saveas(fig_h, figure_name, 'fig')
    print('-dtiff','-r300',[figure_name,'.tif'])
    %print('-dpdf',[fullfile(folder_name,out_file_name),'.pdf'])
    print('-dmeta',[figure_name,'.emf'])
end
%}
end
