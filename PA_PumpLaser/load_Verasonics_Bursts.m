function [PA_data, Bmode_data, burst_info] = load_Verasonics_Bursts(folder) 
%
% function [PA_data, Bmode_data, burst_info] = load_Verasonics_Bursts(folder) 
% 
% The objective of this function is to load data collected using the
% Verasonics system. The data must be stored as bursts (folders names:
% BurstSet###) in the "folder" given as input. 
% 
% INPUTS: 
%   folder - folder containing burst folders (Type: string)
% 
% OUTPUTS: 
%   PA_data    - I/Q data for PA images for all bursts (Type: struct)
%   Bmode_data - I/Q data for Bmode images for all bursts (Type: struct)
%   burst_info - info for the respective burst 
%
% Author: Vinoin Devpaul Vincely (2024)

% Edit Log:
% 02/18/2024 - if/else for direct Acqustion file extraction 
% 07/03/2024 - Included a comprehensive helper function 


current_dir = pwd; 
cd(folder); 

% burst_info = get_Verasonics_BurstsInfo('.'); 
try
    burst_info = get_Verasonics_BurstsInfo('.'); 
catch 
    burst_info = [];
    disp('No burst files present!');
end 

% keyboard
burst_list = ls('./BurstSet*'); 

if isempty(burst_list)
    fprintf('\n Loading acquistion files from %s folder...', folder); 
    acq_list = ls('./*Acq*'); 
    for c = 1:length(acq_list(:,1))
        load_file = load(acq_list(c,:)); 
        PA_data.IData(:,:,c) = load_file.PAIData; 
        PA_data.QData(:,:,c) = load_file.PAQData; 
        Bmode_data.IData(:,:,c) = load_file.BModeIData; 
        Bmode_data.QData(:,:,c) = load_file.BModeQData; 
        clear('load_file'); 
    end 
    fprintf('Complete!\n');
    
else
    fprintf('\nLoading Verasonics Bursts! Loading %d Bursts...\n', length(burst_list(:,1))); 
    
    for b = 1:length(burst_list(:,1))
        cd(burst_list(b,:)); 
        acq_list = ls('./*Acq*'); 
        for c = 1:length(acq_list(:,1))
            load_file = load(acq_list(c,:)); 
            PA_data(b).IData(:,:,c) = load_file.PAIData; 
            PA_data(b).QData(:,:,c) = load_file.PAQData; 
            Bmode_data(b).IData(:,:,c) = load_file.BModeIData; 
            Bmode_data(b).QData(:,:,c) = load_file.BModeQData; 
            clear('load_file'); 
        end 
        
        cd('..'); 
        fprintf(' === Loaded Burst #%d\n', b'); 
        
    end 
    
    fprintf('\nFinished Loading Dataset! \n\n'); 
end 

cd(current_dir);