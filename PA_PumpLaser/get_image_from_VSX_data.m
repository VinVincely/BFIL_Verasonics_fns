function [PA_Bmode] = get_image_from_VSX_data(PA_Bmode, key)
% 
% function [PA_Bmode] = get_image_from_VSX_data(PA_Bmode, key)
% 
% This function retrives images from I/Q data. Data must be placed in a
% structure containing fields "IData" & "QData" (for simplicity use outputs
% from the "load_Verasonics_Bursts" function). 
% 
% INPUTS: 
%   PA_Bmode - I/Q Data for PA/Bmode (Type: struct)
%   key - A key indicating mode of image generation (Type: int)
%       "1" - simple absolute of I/Q data, i.e. img = sqrt(I^2 + Q^2)
%       "2" - log based image construction (NOTE: NEEDS ATTENTION!)
% 
% OUTPUTS:
%   PA_Bmode - Returned input structure with "img" field containing the
%       image generated with I/Q data (Type: struct)
% 
% Author: Vinoin Devpaul Vincely (2024)

% Edit Log: 
%   07/03/2024 - Included a detailed helper funciton

switch key
    
    % A simple image extraction (sqrt(I^2 + Q^2))
    case 1 
        for c = 1:length(PA_Bmode)
            PA_Bmode(c).img = image_from_IQ_simple(PA_Bmode(c).IData, PA_Bmode(c).QData); 
        end
    
    % Log base image extraction
    case 2 
        for c = 1:length(PA_Bmode)
            PA_Bmode(c).img = envelopeDet_logCompression(PA_Bmode(c).IData, PA_Bmode(c).QData); 
        end 
end 