function [exp_para, exp_struct] = experimental_condition(exp_para)
%{
This function create the training protocols of experments.
exp_para: if the input exp_para is [], this function output a default
protocol.
Then the users can modify the parameters in exp_para and run this function
again to pack the parameters into exp_struct.
exp_struct: the standard format of the experimental protocol used in the
simulation function Dx_steady_state_MBON_0301_2023().
%}

%% set default value
if ~isfield(exp_para, 'session')
    % the parameters of each type of sessions
    exp_para.session = ...
        struct('name',   {'test';'imaging';'training'}, ...
               't_CS_plus' ,{  5;        5;       30}, ...
               't_CS_minus',{  5;        5;       30}, ...
               't_ISI'     ,{ 60;      120;      135}, ...
               'n'         ,{  1;        1;        3});
end

%---rest session
if ~isfield(exp_para, 'rest_t')
    % The length of all resting sessions is stored in the rest_t column vector.
    exp_para.rest_t = 3600; 
end
%---session list
if ~isfield(exp_para, 'session_list')
    % The list of all the sessions
    % This function creates the exp_struct by stacking all the sessions by
    % the list.
    exp_para.session_list = ...
        {'imaging';
         'training';
         'imaging';
         'training';
         'imaging';
         'rest';
         'imaging'};
end
%% set session_cell
session_cell = cell(size(exp_para.session,1),1);
for session_i = 1:size(session_cell,1)
    t_length_curr = {exp_para.session(session_i).t_CS_plus; ...
                     exp_para.session(session_i).t_ISI; ...
                     exp_para.session(session_i).t_CS_plus; ...
                     exp_para.session(session_i).t_ISI};
	if strcmp(exp_para.session(session_i).name,'training')
        punishment_curr = {1;0;0;0};
    else
        punishment_curr = {0;0;0;0};
    end
    
    %{imaging_curr = {1;0;2;0};
    if any(strcmp(exp_para.session(session_i).name,{'imaging','extinction'}))
        imaging_curr = {1;0;2;0};
        %specify the matrix index to store the imaging date.
    else
        imaging_curr = {0;0;0;0};
    end
    %}
    
    session_cell_1 = ...
        struct('name'      ,exp_para.session(session_i).name, ...
               'name_sub'  ,{'CS+';'ISI';'CS-';'ISI'}, ...
               't_length'  ,t_length_curr, ...
               'odor'      ,{[1;0];[0;0];[0;1];[0;0]}, ...
               'punishment',punishment_curr, ...
               'imaging_i' ,imaging_curr);
           
    session_cell{session_i,1} = repmat(session_cell_1, exp_para.session(session_i).n,1);
    
end
session_cell = [session_cell; num2cell(...
        struct('name'      ,'rest', ...
               'name_sub'  ,' ', ...
               't_length'  ,num2cell(exp_para.rest_t), ...
               'odor'      ,[0; 0], ...
               'punishment',0, ...
               'imaging_i' ,0))];
%% 
session_name = {exp_para.session.name}';
session_index = zeros(size(exp_para.session_list));
for name_i = 1:length(session_name)
    session_index(strcmp(exp_para.session_list,session_name{name_i})) = name_i;
end
session_index(strcmp(exp_para.session_list,'rest')) = (name_i+1):length(session_cell);

exp_struct = cell2mat(session_cell(session_index));
end