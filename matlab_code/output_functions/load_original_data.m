function [Dx_mean, Dx_SEM] = load_original_data(data_ori_filename, is_simplified_model)
%{
This function read data from data_ori_filename.
% if is_simplified_model == true : 2-module model
% if is_simplified_model == false: 3-module model
%}
data_ori = cat(4, xlsread(data_ori_filename,'ACVvsETA'), ...
                  xlsread(data_ori_filename,'OCTvsBEN'));

%{
Dimension 1: PPL1-v1pedc, PPL1-a'2a2, PPL1-a3, MBON-v1pedc, MBON-a2sc, MBON-a3
Dimension 2: CS+, CS-
Dimension 3: Pre-traing, after 3x training, after 6x training, 1hr memory, 3hr, 24hr
Dimension 4: ACV vs ETA, OCT vs BEN
%}
Dx_mean = permute(cat(3,data_ori(1:3:end,1:6,:,:),data_ori(1:3:end,9:14,:,:)),[1 3 2 4]);
Dx_SEM  = permute(cat(3,data_ori(2:3:end,1:6,:,:),data_ori(2:3:end,9:14,:,:)),[1 3 2 4]);

if is_simplified_model
    % not include PPL1-a'2a2 MBON-a2sc data
    Dx_mean([2 5],:,:,:) = nan;
    Dx_SEM([2 5],:,:,:) = nan;
end
%{
% not include the 24hr data
Dx_mean = Dx_mean(:,:,1:end-1,:);
Dx_SEM = Dx_SEM(:,:,1:end-1,:);
%}
end