clear all
clc
config_path_and_function
% load('Dx_steady_state_nonlinear_3_27-Mar-2023_3modules') %% Fig5i
load('Dx_steady_state_nonlinear_3_27-Mar-2023_2modules') %% Extended Data Fig10m


%%
figure_type_list = {'ACV_EtA','OCT_BEN'};

for figure_type_i= 1:length(figure_type_list)
figure_type = figure_type_list{figure_type_i};
switch figure_type
    case 'ACV_EtA'
        odor_valence = [-3.42;-4.58;-5.60]; % ACV/EtA

    case 'OCT_BEN'
        odor_valence = [2.57;2.24;3.35]; % BEN/OCT
end

% para_mu([12,13,15,18,19],1) = 0; % disconnect MBON-v1 feedbacks; for Extended Data Fig10h
para_mat_cell = reset_parameter_by_valance(para_mu,mat_lu_cell,odor_valence);

shock_strength = 1.5; % shock intensity
para_mat_cell{3}(1) = para_mat_cell{3}(1)*shock_strength;

t_extinction_list = [nan 10 20 30 60 90 120]*60;
rest_t_list = (300:300:3600*3)-135;

output_data_matrix = nan(length(t_extinction_list),length(rest_t_list),3,6);

%%
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
            [Dx_DAN_MBON,data_all] = Dx_steady_state_MBON_revision(para_mat_cell,exp_struct);
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

% figure
% J = customcolormap([0 0.25 0.375 0.5 0.625 0.75 1], {'#8B008B','#C71585','#FA8072','#FFFAFA','#48D1CC','#4682B4','#000000'});
% 
% colormap(J);
% axis_h = zeros(3,6);
% for row_i = 1:3
%     for col_i = 1:6
%         axis_h(row_i,col_i) = subplot(3,6,(row_i-1)*6+col_i); hold on;
%         imagesc(rest_t_list/3600,1:length(t_extinction_list), ...
%                 output_data_matrix(:,:,row_i,col_i),[-40,40])
%     end
% end

switch figure_type
    case 'ACV_EtA'
        output_data_matrix_attractive = output_data_matrix; % ACV/EtA

    case 'OCT_BEN'
        output_data_matrix_repulsive = output_data_matrix; % BEN/OCT
end
end

%%
figure
J = customcolormap([0 0.25 0.375 0.5 0.625 0.75 1], {'#8B008B','#C71585','#FA8072','#FFFAFA','#48D1CC','#4682B4','#000000'});

colormap(J);
subplot(3,2,1)
imagesc(rest_t_list/3600,1:length(t_extinction_list), ...
    output_data_matrix_attractive(:,:,1,6)-Dx_DAN_MBON(6,1,1),[-40,40])
xlabel('Time (hr)','FontSize', 15)
set(gca,'YTickLabel',[]);
title('CS+ induced (Attractive odors)')


subplot(3,2,2)
imagesc(rest_t_list/3600,1:length(t_extinction_list), ...
    output_data_matrix_repulsive(:,:,1,6)-Dx_DAN_MBON(6,1,1),[-40,40])
xlabel('Time (hr)','FontSize', 15)
set(gca,'YTickLabel',[]);
title('CS+ induced (Repulsive odors)')


subplot(3,2,3)
imagesc(rest_t_list/3600,1:length(t_extinction_list), ...
    output_data_matrix_attractive(:,:,2,6)-Dx_DAN_MBON(6,1,1),[-40,40])
xlabel('Time (hr)','FontSize', 15)
set(gca,'YTickLabel',[]);
title('CS– induced (Attractive odors)')


subplot(3,2,4)
imagesc(rest_t_list/3600,1:length(t_extinction_list), ...
    output_data_matrix_repulsive(:,:,2,6)-Dx_DAN_MBON(6,1,1),[-40,40])
xlabel('Time (hr)','FontSize', 15)
set(gca,'YTickLabel',[]);
title('CS– induced (Repulsive odors)')


subplot(3,2,5)
imagesc(rest_t_list/3600,1:length(t_extinction_list), ...
    output_data_matrix_attractive(:,:,3,6),[-40,40])
xlabel('Time (hr)','FontSize', 15)
set(gca,'YTickLabel',[]);
title('CS+ vs. CS– induced (Attractive odors)')


subplot(3,2,6)
imagesc(rest_t_list/3600,1:length(t_extinction_list), ...
    output_data_matrix_repulsive(:,:,3,6),[-40,40])
xlabel('Time (hr)','FontSize', 15)
set(gca,'YTickLabel',[]);
title('CS+ vs. CS– induced (Repulsive odors)')

