function [nAcq_PA_Bmode] = get_nAcq_VSX_Bursts(PA_Bmode, n)
% 
% function [nAcq_PA_Bmode] = get_nAcq_VSX_Bursts(PA_Bmode, n)
% 
% The objective of this function is to get the "n-th" acquisiton from the
% VSX Bursts structure.
% 
% INPUTS: 
%   PA_Bmode - I/Q or image data for PA/Bmode (Type: struct)
%   n        - Acquistion numbers to be retrived (Type: vector/scalar int)
% 
% OUTPUT:
%   PA_Bmode - I/Q or image data for PA/Bmode with desired acquistions 
%       extracted (Type: struct)
%
% Author: Vinoin Devpaul Vincely (2024)

% Edit History 
%   04/01/2024 - included if/else for RF data input
%   07/03/2024 - Included a detailed helper funciton

% ==== User notification 
if isfield(PA_Bmode, 'data')
    fprintf('\tReturning RF data of Acq. #%d!\n', n); 
else 
    if ~isfield(PA_Bmode, 'img')
        fprintf('\tReturning I & Q data of Acq. #%d!\n', n); 
    end 
end

for c = 1:length(PA_Bmode)
    x = 1;
    for acq = n
        
        if isfield(PA_Bmode, 'data')
            nAcq_PA_Bmode(c).data(:,:,x) = PA_Bmode(c).data(:,:,acq);  
        else 
            if acq > size(PA_Bmode(c).img,3)
    %             keyboard
                continue; 
            else 
                if isfield(PA_Bmode, 'img')
                    nAcq_PA_Bmode(c).img(:,:,x) = PA_Bmode(c).img(:,:,acq); %keyboard
                else 
                    nAcq_PA_Bmode(c).IData(:,:,x) = PA_Bmode(c).IData(:,:,acq);
                    nAcq_PA_Bmode(c).QData(:,:,x) = PA_Bmode(c).QData(:,:,acq);
                end 
            end            
        end 
        x = x + 1; 
    end 
end 
