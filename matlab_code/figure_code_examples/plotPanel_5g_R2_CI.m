clear all
clc
config_path_and_function
load('Dx_steady_state_nonlinear_3_27-Mar-2023_3modules') % Fig5g
% load('Dx_steady_state_nonlinear_3_27-Mar-2023_2modules') %% Extended Data Fig10l


MemoryTrace_CSp = cell(6,3);
MemoryTrace_CSm = cell(6,3);
MemoryTrace_delta = cell(6,3);

for i = 1:length(para_rand)

    odor_valence = 0;
    para_mat_cell = reset_parameter_by_valance(para_rand(:,i),mat_lu_cell,odor_valence);

    %-----loading model fitted parameters
    [exp_para, exp_struct] = experimental_condition([]);
    [exp_para, exp_struct] = experimental_condition(exp_para);

    exp_para.session_list = ...
        {'imaging';
        'training';
        'rest';
        'imaging'};

    shock_strength = 1; % shock intensity
    para_mat_cell{3}(1) = para_mat_cell{3}(1)*shock_strength;

    exp_para.session(3).n = 10; % training cycle
    exp_para.session(3).t_ISI = 135; % interval

    ISI_list = [60,135,360,600,900];
    rest_t_list = [300,10800,86400];

    [exp_para, exp_struct] = experimental_condition(exp_para);


    Dx_DAN_MBON_Matrix_CSp = cell(6,1);
    Dx_DAN_MBON_Matrix_CSm = cell(6,1);
    Dx_DAN_MBON_Matrix_delta = cell(6,1);

    indx = 1;
    for ISI_n = ISI_list

        exp_para.session(3).t_ISI = ISI_n; % interval
        exp_struct(4*(exp_para.session(3).n+1)).t_length = 0;

        indxx = 1;

        for rest_t = rest_t_list

            exp_para.rest_t = rest_t;
            [exp_para, exp_struct] = experimental_condition(exp_para);
            exp_struct(4).t_length = 300;
            [Dx_DAN_MBON{i},data_all] = Dx_steady_state_MBON_revision(para_mat_cell,exp_struct);

            for neuron_inx = 1:6
                for ind = 1:size(Dx_DAN_MBON{i},3)
                    Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON{i}(neuron_inx,1,ind);
                    Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON{i}(neuron_inx,2,ind);
                    Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON{i}(neuron_inx,1,ind) - Dx_DAN_MBON{i}(neuron_inx,2,ind);
                end
            end
            indxx = indxx+1;
        end
        for neuron_inx = 1:6
            MemoryTrace_CSp{neuron_inx,indx}(i,1:2) = Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(1,1:2);
            MemoryTrace_CSm{neuron_inx,indx}(i,1:2) = Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(1,1:2);
            MemoryTrace_delta{neuron_inx,indx}(i,1:2) = Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(1,1:2);

            MemoryTrace_CSp{neuron_inx,indx}(i,2:indxx) = Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(:,size(Dx_DAN_MBON{i},3))';
            MemoryTrace_CSm{neuron_inx,indx}(i,2:indxx) = Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(:,size(Dx_DAN_MBON{i},3))';
            MemoryTrace_delta{neuron_inx,indx}(i,2:indxx) = Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(:,size(Dx_DAN_MBON{i},3))';
        end
        indx = indx+1;
    end
end



%%
for i = 1:5
    for neuron_inx = 1:6
        MemoryTrace_delta_5min{neuron_inx}(1:10000,i)= MemoryTrace_delta{neuron_inx,i}(:,2);
        MemoryTrace_delta_3hr{neuron_inx}(1:10000,i)= MemoryTrace_delta{neuron_inx,i}(:,3);
        MemoryTrace_delta_24hr{neuron_inx}(1:10000,i)= MemoryTrace_delta{neuron_inx,i}(:,4);

        MemoryTrace_CSp_5min{neuron_inx}(1:10000,i)= MemoryTrace_CSp{neuron_inx,i}(:,2);
        MemoryTrace_CSp_3hr{neuron_inx}(1:10000,i)= MemoryTrace_CSp{neuron_inx,i}(:,3);
        MemoryTrace_CSp_24hr{neuron_inx}(1:10000,i)= MemoryTrace_CSp{neuron_inx,i}(:,4);

        MemoryTrace_CSm_5min{neuron_inx}(1:10000,i)= MemoryTrace_CSm{neuron_inx,i}(:,2);
        MemoryTrace_CSm_3hr{neuron_inx}(1:10000,i)= MemoryTrace_CSm{neuron_inx,i}(:,3);
        MemoryTrace_CSm_24hr{neuron_inx}(1:10000,i)= MemoryTrace_CSm{neuron_inx,i}(:,4);

        %------------------------------------------------------------------
        A = [];
        A = MemoryTrace_delta_5min{neuron_inx}(1:10000,i);
        A = A(~isnan(A));
        A = A(~isinf(A));
        MemoryTrace_delta_5min_avg(neuron_inx,i) = median(A);

        A = sort(A);
        p = 0.68;
        MemoryTrace_delta_5min_CI{neuron_inx,1}(1,i) = A(round(length(A)*(1+p)/2))-median(A);
        MemoryTrace_delta_5min_CI{neuron_inx,1}(2,i) = median(A) - A(round(length(A)*(1-p)/2));

        AA = [];
        AA = MemoryTrace_delta_3hr{neuron_inx}(1:10000,i);
        AA = AA(~isnan(AA));
        AA = AA(~isinf(AA)); 
        MemoryTrace_delta_3hr_avg(neuron_inx,i) = median(AA);
        

        AA = sort(AA);
        p = 0.68;
        MemoryTrace_delta_3hr_CI{neuron_inx,1}(1,i) = AA(round(length(AA)*(1+p)/2))-median(AA);
        MemoryTrace_delta_3hr_CI{neuron_inx,1}(2,i) = median(AA) - AA(round(length(AA)*(1-p)/2));

        AAA = [];
        AAA = MemoryTrace_delta_24hr{neuron_inx}(1:10000,i);
        AAA = AAA(~isnan(AAA));
        AAA = AAA(~isinf(AAA)); 
        MemoryTrace_delta_24hr_avg(neuron_inx,i) = median(AAA);
       

        AAA = sort(AAA);
        p = 0.68;
        MemoryTrace_delta_24hr_CI{neuron_inx,1}(1,i) = AAA(round(length(AAA)*(1+p)/2))-median(AAA);
        MemoryTrace_delta_24hr_CI{neuron_inx,1}(2,i) = median(AAA) - AAA(round(length(AAA)*(1-p)/2));

        %------------------------------------------------------------------
        A = [];
        A = MemoryTrace_CSp_5min{neuron_inx}(1:10000,i);
        A = A(~isnan(A));
        A = A(~isinf(A));
        MemoryTrace_CSp_5min_avg(neuron_inx,i) = median(A);

        A = sort(A);
        p = 0.68;
        MemoryTrace_CSp_5min_CI{neuron_inx,1}(1,i) = A(round(length(A)*(1+p)/2))-median(A);
        MemoryTrace_CSp_5min_CI{neuron_inx,1}(2,i) = median(A) - A(round(length(A)*(1-p)/2));

        AA = [];
        AA = MemoryTrace_CSp_3hr{neuron_inx}(1:10000,i);
        AA = AA(~isnan(AA));
        AA = AA(~isinf(AA));
        MemoryTrace_CSp_3hr_avg(neuron_inx,i) = median(AA);

        AA = sort(AA);
        p = 0.68;
        MemoryTrace_CSp_3hr_CI{neuron_inx,1}(1,i) = AA(round(length(AA)*(1+p)/2))-median(AA);
        MemoryTrace_CSp_3hr_CI{neuron_inx,1}(2,i) = median(AA) - AA(round(length(AA)*(1-p)/2));

        AAA = [];
        AAA = MemoryTrace_CSp_24hr{neuron_inx}(1:10000,i);
        AAA = AAA(~isnan(AAA));
        AAA = AAA(~isinf(AAA));
        MemoryTrace_CSp_24hr_avg(neuron_inx,i) = median(AAA);

        AAA = sort(AAA);
        p = 0.68;
        MemoryTrace_CSp_24hr_CI{neuron_inx,1}(1,i) = AAA(round(length(AAA)*(1+p)/2))-median(AAA);
        MemoryTrace_CSp_24hr_CI{neuron_inx,1}(2,i) = median(AAA) - AAA(round(length(AAA)*(1-p)/2));

        %------------------------------------------------------------------
        A = [];
        A = MemoryTrace_CSm_5min{neuron_inx}(1:10000,i);
        A = A(~isnan(A));
        A = A(~isinf(A));
        MemoryTrace_CSm_5min_avg(neuron_inx,i) = median(A);
        A = sort(A);
        p = 0.68;
        MemoryTrace_CSm_5min_CI{neuron_inx,1}(1,i) = A(round(length(A)*(1+p)/2))-median(A);
        MemoryTrace_CSm_5min_CI{neuron_inx,1}(2,i) = median(A) - A(round(length(A)*(1-p)/2));

        AA = [];
        AA = MemoryTrace_CSm_3hr{neuron_inx}(1:10000,i);
        AA = AA(~isnan(AA));
        AA = AA(~isinf(AA));
        MemoryTrace_CSm_3hr_avg(neuron_inx,i) = median(AA);
        AA = sort(AA);
        p = 0.68;
        MemoryTrace_CSm_3hr_CI{neuron_inx,1}(1,i) = AA(round(length(AA)*(1+p)/2))-median(AA);
        MemoryTrace_CSm_3hr_CI{neuron_inx,1}(2,i) = median(AA) - AA(round(length(AA)*(1-p)/2));

        AAA = [];
        AAA = MemoryTrace_CSm_24hr{neuron_inx}(1:10000,i);
        AAA = AAA(~isnan(AAA));
        AAA = AAA(~isinf(AAA));
        MemoryTrace_CSm_24hr_avg(neuron_inx,i) = median(AAA);
        AAA = sort(AAA);
        p = 0.68;
        MemoryTrace_CSm_24hr_CI{neuron_inx,1}(1,i) = AAA(round(length(AAA)*(1+p)/2))-median(AAA);
        MemoryTrace_CSm_24hr_CI{neuron_inx,1}(2,i) = median(AAA) - AAA(round(length(AAA)*(1-p)/2));

    end
end

%% Plot MBONs
for MBON_inx = 1:3
    switch MBON_inx
        case 1
            MBON_v1_diff = [MemoryTrace_delta_5min_avg(4,:), 0, 0, MemoryTrace_delta_3hr_avg(4,:), 0, 0, MemoryTrace_delta_24hr_avg(4,:)];
            MBON_v1_diff(2:3,1:5) = MemoryTrace_delta_5min_CI{4,1};
            MBON_v1_diff(2:3,8:12) = MemoryTrace_delta_3hr_CI{4,1};
            MBON_v1_diff(2:3,15:19) = MemoryTrace_delta_24hr_CI{4,1};

            data = MBON_v1_diff(1,:);
            errhigh = MBON_v1_diff(2,:);
            errlow  = MBON_v1_diff(3,:);
            MBON_name = 'MBON-\gamma1pedc';
        case 2
            MBON_a2_diff = [MemoryTrace_delta_5min_avg(5,:), 0, 0, MemoryTrace_delta_3hr_avg(5,:), 0, 0, MemoryTrace_delta_24hr_avg(5,:)];
            MBON_a2_diff(2:3,1:5) = MemoryTrace_delta_5min_CI{5,1};
            MBON_a2_diff(2:3,8:12) = MemoryTrace_delta_3hr_CI{5,1};
            MBON_a2_diff(2:3,15:19) = MemoryTrace_delta_24hr_CI{5,1};

            data = MBON_a2_diff(1,:);
            errhigh = MBON_a2_diff(2,:);
            errlow  = MBON_a2_diff(3,:);
            MBON_name = 'MBON-\alpha2';
        case 3
            MBON_a3_diff = [MemoryTrace_delta_5min_avg(6,:), 0, 0, MemoryTrace_delta_3hr_avg(6,:), 0, 0, MemoryTrace_delta_24hr_avg(6,:)];
            MBON_a3_diff(2:3,1:5) = MemoryTrace_delta_5min_CI{6,1};
            MBON_a3_diff(2:3,8:12) = MemoryTrace_delta_3hr_CI{6,1};
            MBON_a3_diff(2:3,15:19) = MemoryTrace_delta_24hr_CI{6,1};

            data = MBON_a3_diff(1,:);
            errhigh = MBON_a3_diff(2,:);
            errlow  = MBON_a3_diff(3,:);
            MBON_name = 'MBON-\alpha3';
    end

    subplot(1,3,MBON_inx)
    X = 1:19;
    bar(X,data)
    set(gca,'XTickLabel',[]);
    title(MBON_name,'FontSize', 15)
    box off
    hold on


    er = errorbar(X,data,errlow,errhigh);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    ylabel('CS+ vs. CS– bias','FontSize', 15)

    text(1,-46,'5-min','FontSize', 15)
    text(9,-46,'3-hr','FontSize', 15)
    text(15,-46,'24-hr','FontSize', 15)


    ylim([-45 15])
    hold off

end