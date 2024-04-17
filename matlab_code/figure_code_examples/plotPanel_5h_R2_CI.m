 clear all
clc
config_path_and_function
load('Dx_steady_state_nonlinear_3_27-Mar-2023_3modules')

MemoryTrace_CSp = cell(6,5);
MemoryTrace_CSm = cell(6,5);
MemoryTrace_delta = cell(6,5);
ISI_list = [60,135,360,600,900];

indx = 1;
for ISI_n = ISI_list
    for i = 1:length(para_rand)

        odor_valence = [-3.42;-4.58;-5.60]; % ACV/EtA
        % odor_valence = [2.57;2.24;3.35]; % BEN/OCT
        % odor_valence = 0;

        shock_strength = 1.5; % shock intensity
        para_mat_cell = reset_parameter_by_valance(para_rand(:,i),mat_lu_cell,odor_valence);
        para_mat_cell{3}(1) = para_mat_cell{3}(1)*shock_strength;


        %-----loading model fitted parameters
        [exp_para, exp_struct] = experimental_condition([]);
        [exp_para, exp_struct] = experimental_condition(exp_para);

        exp_para.session_list = ...
            {'imaging';
            'training';
            'imaging';
            'rest';
            'imaging'};

        exp_para.session(3).n = 6; % training cycle
        [exp_para, exp_struct] = experimental_condition(exp_para);

        Dx_DAN_MBON_Matrix_CSp = cell(6,1);
        Dx_DAN_MBON_Matrix_CSm = cell(6,1);
        Dx_DAN_MBON_Matrix_delta = cell(6,1);

        indxx = 1;
        exp_para.rest_t = 86400-300-250;
        exp_para.session(3).t_ISI = ISI_n; % interval
        [exp_para, exp_struct] = experimental_condition(exp_para);
        exp_struct(4).t_length = 300;
        exp_struct(4*(exp_para.session(3).n+1)).t_length = 300;

        [Dx_DAN_MBON{i},data_all] = Dx_steady_state_MBON_revision(para_mat_cell,exp_struct);

        for neuron_inx = 1:6
            for ind = 1:size(Dx_DAN_MBON{i},3)
                Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON{i}(neuron_inx,1,ind);
                Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON{i}(neuron_inx,2,ind);
                Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON{i}(neuron_inx,1,ind) - Dx_DAN_MBON{i}(neuron_inx,2,ind);
            end
        end
        indxx = indxx+1;
    
    for neuron_inx = 1:6
        MemoryTrace_CSp{neuron_inx,indx}(i,1:2) = Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(1,1:2);
        MemoryTrace_CSm{neuron_inx,indx}(i,1:2) = Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(1,1:2);
        MemoryTrace_delta{neuron_inx,indx}(i,1:2) = Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(1,1:2);
        
        MemoryTrace_CSp{neuron_inx,indx}(i,3:indxx+1) = Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(:,size(Dx_DAN_MBON{i},3))';
        MemoryTrace_CSm{neuron_inx,indx}(i,3:indxx+1) = Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(:,size(Dx_DAN_MBON{i},3))';
        MemoryTrace_delta{neuron_inx,indx}(i,3:indxx+1) = Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(:,size(Dx_DAN_MBON{i},3))';
    end

    end
    indx = indx+1;
end



%% Median and CI
for i = 1:5
    for neuron_inx = 1:6
        MemoryTrace_CSp_avg{neuron_inx,i} = mean(MemoryTrace_CSp{neuron_inx,i}(:,:));
        MemoryTrace_CSp_sem{neuron_inx,i} = std((MemoryTrace_CSp{neuron_inx,i}(:,:)))/sqrt(14);
        MemoryTrace_CSp_sdt{neuron_inx,i} = std((MemoryTrace_CSp{neuron_inx,i}(:,:)));

        MemoryTrace_CSm_avg{neuron_inx,i} = mean(MemoryTrace_CSm{neuron_inx,i}(:,:));
        MemoryTrace_CSm_sem{neuron_inx,i} = std((MemoryTrace_CSm{neuron_inx,i}(:,:)))/sqrt(14);
        MemoryTrace_CSm_sdt{neuron_inx,i} = std((MemoryTrace_CSm{neuron_inx,i}(:,:)));

        MemoryTrace_delta_5min{neuron_inx}(1:10000,i)= MemoryTrace_delta{neuron_inx,i}(:,2);
        A = [];
        A = MemoryTrace_delta_5min{neuron_inx}(1:10000,i);
        A = A(~isnan(A));
        A = A(~isinf(A));
        MemoryTrace_delta_5min_avg(neuron_inx,i) = median(A);

        A = sort(A);
        p = 0.68;
        MemoryTrace_delta_5min_CI{neuron_inx,1}(1,i) = A(round(length(A)*(1+p)/2))-median(A);
        MemoryTrace_delta_5min_CI{neuron_inx,1}(2,i) = median(A) - A(round(length(A)*(1-p)/2));

        MemoryTrace_delta_24hr{neuron_inx}(1:10000,i)= MemoryTrace_delta{neuron_inx,i}(:,3);
        A = [];
        A = MemoryTrace_delta_24hr{neuron_inx}(1:10000,i);
        A = A(~isnan(A));
        A = A(~isinf(A));
        MemoryTrace_delta_24hr_avg(neuron_inx,i) = median(A);

        A = sort(A);
        p = 0.68;
        MemoryTrace_delta_24hr_CI{neuron_inx,1}(1,i) = A(round(length(A)*(1+p)/2))-median(A);
        MemoryTrace_delta_24hr_CI{neuron_inx,1}(2,i) = median(A) - A(round(length(A)*(1-p)/2));
    end
end

%% plot MBON-a3
MBON_a3_diff = [MemoryTrace_delta_5min_avg(6,[1,3,5]), 0, 0, MemoryTrace_delta_24hr_avg(6,[1,3,5])];
MBON_a3_diff(2:3,1:3) = MemoryTrace_delta_5min_CI{6,1}(:,[1,3,5]);
MBON_a3_diff(2:3,6:8) = MemoryTrace_delta_24hr_CI{6,1}(:,[1,3,5]);

data = MBON_a3_diff(1,:);
errhigh = MBON_a3_diff(2,:);
errlow  = MBON_a3_diff(3,:);
MBON_name = 'MBON-\alpha3';

X = 1:8;
bar(X,data)
set(gca,'XTickLabel',[]);
title(MBON_name,'FontSize', 15)
box off
hold on


er = errorbar(X,data,errlow,errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylabel('CS+ vs. CS– bias','FontSize', 15)

text(2,-46,'5-min','FontSize', 15)
text(7,-46,'24-hr','FontSize', 15)
text(0.8,5,'60s-ISI','FontSize', 15)
text(1.8,5,'360s-ISI','FontSize', 15)
text(2.8,5,'900s-ISI','FontSize', 15)
text(5.8,5,'60s-ISI','FontSize', 15)
text(6.8,5,'360s-ISI','FontSize', 15)
text(7.8,5,'900s-ISI','FontSize', 15)


ylim([-45 10])
hold off
