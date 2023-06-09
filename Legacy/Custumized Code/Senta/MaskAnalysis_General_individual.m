%% Covid-19 Slide Analysis for Janielle
% 08/10/2021
% Peiyu Wang

% For Janielle only, the mask can be read by filefolder/ROI_Mask/ControlMask.tif;

% Data only has one z plane. which is pretty good! 

% Let's do this. 

close all; clear all;
addpath("Functions");

DataFolder = "D:\Scotts Lab\Collaborations\For Senta\Processed Pics";
%% Tunable Variables 
mode_iter = 10;
thresh_val = 3; 

%% Predifining Variables

condition = [];
islet_No = [];
individual = [];
G_mode_mask = [];
S_mode_mask = [];

G_mode_full = [];
S_mode_full = [];

cond_folder = dir(DataFolder);
for cond_idx = 3:numel(cond_folder)
    ind_folder = dir(fullfile(DataFolder,cond_folder(cond_idx).name));
    for ind_idx = 3: numel(ind_folder)
        samp_folder = dir(fullfile(ind_folder(ind_idx).folder,ind_folder(ind_idx).name));
        for samp_idx = 3: numel(samp_folder)
            data_folder = fullfile(samp_folder(samp_idx).folder,samp_folder(samp_idx).name);
            
            condition = cat(1,condition,string(cond_folder(cond_idx).name));
            individual = cat(1,individual,string(ind_folder(ind_idx).name));
            islet_No = cat(1,islet_No,string(samp_folder(samp_idx).name));
                      
            imageFile = dir(fullfile(data_folder,'*.tif'));
            int = imread(fullfile(imageFile(9).folder,imageFile(9).name));
            G = standardPhase(imread(fullfile(imageFile(11).folder,imageFile(11).name)));
            S = standardPhase(imread(fullfile(imageFile(12).folder,imageFile(12).name)));
            mask_img = imread(fullfile(data_folder,"ROI_Mask","ControlMask.tif"));  %% Change the name here to the name of the mask
            
            current_struct = struct('int',int,'G',G,'S',S);
            current_struct = threshstruct(current_struct,thresh_val);
            [G_f,S_f] = findModePhasor(current_struct,mode_iter);
            
            G_mode_full = cat(1,G_mode_full,G_f);
            S_mode_full = cat(1,S_mode_full,S_f);
            
            mask_struct = maskphasorstruct(current_struct,mask_img);
            mask_struct = threshstruct(mask_struct,thresh_val);
            [G_m,S_m] = findModePhasor(mask_struct,mode_iter);
            
            G_mode_mask = cat(1,G_mode_mask,G_m);
            S_mode_mask = cat(1,S_mode_mask,S_m);
        end
    end
end

%% analysis for individuals.
condition_names = unique(condition);
individual_name = unique(individual);

G_mode_sum_mask = [];
S_mode_sum_mask = [];

G_mode_sum_full = [];
S_mode_sum_full = [];

condition_sum = [];
individual_sum = [];

for i = 1: numel(individual_name)
    index = find(individual == individual_name(i));
    G_current = G_mode_mask(index);
    S_current = G_mode_mask(index);
    G_mode_sum_mask = cat(1,G_mode_sum_mask,mean(G_current(:)));
    S_mode_sum_mask = cat(1,S_mode_sum_mask,mean(S_current(:)));
    
    G_current = G_mode_full(index);
    S_current = S_mode_full(index);
    G_mode_sum_full = cat(1,G_mode_sum_full,mean(G_current(:)));
    S_mode_sum_full = cat(1,S_mode_sum_full,mean(S_current(:)));
    
    
    
    condition_sum =  cat(1,condition_sum,condition(index(1)));
    individual_sum =  cat(1,individual_sum, individual_name(i));
end
%% Line Extention Metabolism for Free Precentage Calculation
[mask_Mode_LEXT,mask_G_int,mask_S_int,mask_tao] = lineExtensionMetabolism(G_mode_sum_mask,S_mode_sum_mask);
[full_Mode_LEXT,full_G_int,full_S_int,full_tao] = lineExtensionMetabolism(G_mode_sum_full,S_mode_sum_full);

%%

DataTable=table(condition_sum,individual_sum,mask_Mode_LEXT,full_Mode_LEXT,G_mode_sum_mask,S_mode_sum_mask,G_mode_sum_full,S_mode_sum_full...
    ,mask_G_int,mask_S_int,mask_tao,full_G_int,full_S_int,full_tao);
filename = 'individual.xlsx';
writetable(DataTable,filename,'Sheet',1)

%% Functions
function new_struct = maskphasorstruct(org_struct, map)

int = org_struct.int; G = org_struct.G; S = org_struct.S;
int(map == 0) = 0; G(map == 0) = 0; S(map == 0) = 0;

new_struct = struct('int',int,'G',G,'S',S);
end



%%
%% Function: findModePhasor
% Peiyu Wang
% 11/4/2021

% Finding the phasor mode. 
% This is done by finding a weighted mode on G and S respectively. 
% The number of pixel in the phasor plot taken into consideration is
% determined by mode_iter. 

function [G_mode,S_mode] = findModePhasor(org_ref, mode_iter)


map_res = 1024; 
phasor_his = zeros(map_res,map_res);

for i = 1:size(org_ref.int,1)
    for j = 1:size(org_ref.int,2)
        G_index = floor((org_ref.G(i,j)-1.526e-05)*map_res/2)+map_res/2+1; %function floor is doing the binning for you. 
        S_index = map_res/2 - floor((org_ref.S(i,j)-1.526e-05)*map_res/2);
        if G_index < 1; G_index = 1; end
        if S_index < 1; S_index = 1; end
        if G_index > map_res; G_index = map_res; end
        if S_index > map_res; S_index = map_res; end
        if G_index && S_index; phasor_his(S_index,G_index) = phasor_his(S_index,G_index)+1;end
    end
end

% because the pixel value at (0,0) is too high, we change that to 0;
[max_val,max_Idx] = max(phasor_his(:));
phasor_his(max_Idx) = 0;


[phasor_his_sorted,phasor_his_index] = sort(phasor_his(:));

G_mode = 0;
S_mode = 0;
total_num = 0;
for i = 1: mode_iter
    G_cur_index = ceil(phasor_his_index(end-i+1)/map_res);
    S_cur_index = rem(phasor_his_index(end-i+1),map_res);
    G_mode = G_mode + phasor_his_sorted(end-i+1) * (-1+G_cur_index/map_res*2 - 1/map_res*2);
    S_mode = S_mode + phasor_his_sorted(end-i+1) * (1-S_cur_index/map_res*2);
    total_num = total_num + phasor_his_sorted(end-i+1);
end

G_mode = G_mode/total_num;
S_mode = S_mode/total_num;

end

%% Functions
function thresh_struct = threshstruct(org_struct,thresh_val)
thresh_struct = org_struct;
thresh_struct.int(org_struct.int<thresh_val) = 0;
thresh_struct.G(org_struct.int<thresh_val) = 0;
thresh_struct.S(org_struct.int<thresh_val) = 0;

end
