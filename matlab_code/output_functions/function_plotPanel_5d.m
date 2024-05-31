function function_plotPanel_5d(file_parameter,figure_name,figure_type)


%load('/Users/flyspike/Desktop/Fly leanring paper/fly model with CI/version2 - a2 a3 equal tau/figures/Dx_steady_state_nonlinear_3_18-Aug-2021')
%load('figures/Dx_steady_state_nonlinear_3_18-Aug-2021')
load(file_parameter)
addpath('customcolormap\')

%-----loading model fitted parameters
%para_mat_cell = parameter_vec2mat(x_lu,mat_lu_cell);



% create the defualt experimental condition listexp_para
[exp_para, exp_struct] = experimental_condition([]);

exp_para.session_list = ...
    {'imaging';
    'training';
    'rest';
    'imaging'};

switch figure_type
    case '5d'
        exp_para.session(3).t_ISI = 135; % interval
        exp_para.session(3).n = 1; % training cycle % Figure 5d
        odor_valence_list = -4:1:4;
        rest_t_list = (300:300:3600)-135;
    case 'S9c'
        exp_para.session(3).t_ISI = 135; % interval
        exp_para.session(3).n = 6; % training cycle % S-figure 8c
        odor_valence_list = -8:1:8;
        rest_t_list = (300:300:86400)-135;
    case '2train360'
        exp_para.session(3).t_ISI = 360; % interval
        exp_para.session(3).n = 2; 
        odor_valence_list = -8:1:8;
        rest_t_list = (300:300:86400)-135;
    case '2train3600'
        exp_para.session(3).t_ISI = 3600; % interval
        exp_para.session(3).n = 2; 
        odor_valence_list = -8:1:8;
        rest_t_list = (300:300:86400)-135;
    otherwise
        error('figure_type should be 5d, S9c, 2train360 or 2train3600.')
end
[exp_para, exp_struct] = experimental_condition(exp_para);

%%
Dx_DAN_MBON_Matrix_CSp = cell(6,1);
Dx_DAN_MBON_Matrix_CSm = cell(6,1);
Dx_DAN_MBON_Matrix_delta = cell(6,1);
MemoryTrace_CSp = cell(6,1);
MemoryTrace_CSm = cell(6,1);
MemoryTrace_delta = cell(6,1);

indx = 1;

for odor_valence = odor_valence_list
    %-4:1:4 % Figure 5d
    % for odor_valence = -8:1:8 % S-figure 8c
    %{
    NoValence = W_1*W_KDKM_KC_0';
    NoValence(1:3,:) = 1 * odor_valence*(2/(1+exp(-KC_adapt*5)));
    W_KDKM_KC_0 = (inv(W_1)*NoValence)';
    %}
    
    %{
    %% ----------ACV-ETA
    KC_adapt = 0.05;
    Firing_Temp(4,:) = 31.6*(2/(1+exp(-KC_adapt*5)));
    Firing_Temp(5,:) = 5.8*(2/(1+exp(-KC_adapt*5)));
    Firing_Temp(6,:) = 17.8*(2/(1+exp(-KC_adapt*5)));
    %}
    para_mat_cell = reset_parameter_by_valance(para_mu,mat_lu_cell,odor_valence);
    indxx = 1;
    
    for rest_t = rest_t_list
        exp_para.rest_t = rest_t;
        [exp_para, exp_struct] = experimental_condition(exp_para);
        exp_struct(4).t_length = 300;
        
        [Dx_DAN_MBON,data_all] = Dx_steady_state_MBON_0301_2023(para_mat_cell,exp_struct);
        for neuron_inx = 1:6
            for ind = 1:size(Dx_DAN_MBON,3)
                Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON(neuron_inx,1,ind);
                Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON(neuron_inx,2,ind);
                Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON(neuron_inx,1,ind) - Dx_DAN_MBON(neuron_inx,2,ind);
            end
        end
        indxx = indxx+1;
    end
    for neuron_inx = 1:6
        MemoryTrace_CSp{neuron_inx,1}(indx,1:2) = Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(1,1:2);
        MemoryTrace_CSm{neuron_inx,1}(indx,1:2) = Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(1,1:2);
        MemoryTrace_delta{neuron_inx,1}(indx,1:2) = Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(1,1:2);
        
        MemoryTrace_CSp{neuron_inx,1}(indx,2:indxx) = Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(:,size(Dx_DAN_MBON,3))';
        MemoryTrace_CSm{neuron_inx,1}(indx,2:indxx) = Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(:,size(Dx_DAN_MBON,3))';
        MemoryTrace_delta{neuron_inx,1}(indx,2:indxx) = Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(:,size(Dx_DAN_MBON,3))';
    end
    indx = indx+1;
end

%%
paper_size=[6 11];

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
axis_h = zeros(6,3);
for neuron_inx = 1:3
    %------------------------------------------------------------------------
    axis_h(1,neuron_inx) = subplot(7,3,neuron_inx);
    imagesc((rest_t_list+135)/3600,odor_valence_list, ...
        MemoryTrace_CSp{neuron_inx,1}(:,2:end)-MemoryTrace_CSp{neuron_inx,1}(:,1),[-40,40])
    
    axis_h(2,neuron_inx) = subplot(7,3,3+neuron_inx);
    imagesc((rest_t_list+135)/3600,odor_valence_list, ...
        MemoryTrace_CSm{neuron_inx,1}(:,2:end)-MemoryTrace_CSm{neuron_inx,1}(:,1),[-40,40])
    
    axis_h(3,neuron_inx) = subplot(7,3,neuron_inx+6);
    imagesc((rest_t_list+135)/3600,odor_valence_list, ...
        MemoryTrace_delta{neuron_inx,1}(:,2:end),[-40,40])
end

for neuron_inx = 4:6
    %------------------------------------------------------------------------
    axis_h(4,neuron_inx-3) = subplot(7,3,neuron_inx+9);
    imagesc((rest_t_list+135)/3600,odor_valence_list, ...
        MemoryTrace_CSp{neuron_inx,1}(:,2:end)-MemoryTrace_CSp{neuron_inx,1}(:,1),[-40 40])
    
    axis_h(5,neuron_inx-3) = subplot(7,3,neuron_inx+12);
    imagesc((rest_t_list+135)/3600,odor_valence_list, ...
        MemoryTrace_CSm{neuron_inx,1}(:,2:end)-MemoryTrace_CSm{neuron_inx,1}(:,1),[-40 40])
        
    axis_h(6,neuron_inx-3) = subplot(7,3,neuron_inx+15);
    imagesc((rest_t_list+135)/3600,odor_valence_list, ...
        MemoryTrace_delta{neuron_inx,1}(:,2:end),[-40 40])
end

set(axis_h,'box','off','YDir','normal');%, ...
    %'YTick',1:2:9,'YTickLabel',-4:2:4, ...
    %'XTick',[1 6.5 12],'XTickLabel',0:0.5:1)

subplot(7,3,1);title('DAN-v1')
subplot(7,3,2);title('DAN-a2')
subplot(7,3,3);title('DAN-a3')

subplot(7,3,13);title('MBON-v1')
subplot(7,3,14);title('MBON-a2')
subplot(7,3,15);title('MBON-a3')

%% save figures
%{
saveas(fig_h, figure_name, 'fig')
print('-dtiff','-r300',[figure_name,'.tif'])
%print('-dpdf',[fullfile(folder_name,out_file_name),'.pdf'])
print('-dmeta',[figure_name,'.emf'])
%}

end
