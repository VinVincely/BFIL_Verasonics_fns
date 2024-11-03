function [hp] = visualize_VSX_sPAI_data(PA, BmodeImg, meanFlnc, lambdas, nPlt)
% 
% function [hp] = visualize_VSX_sPAI_data(PA, BmodeImg, meanFlnc, lambdas, nPlt)
% 
% The objective of this function is to visualize the acquired PA/US data. 
% This function can be used to visualize multiple bursts simultaneously. 
% The variables used here were written with the assumption that each burst
% was acquired at a different wavelength. 
% 
% INPUT: 
%   PA - cell of photoacoustic images
%   BmodeImg - cell of B-mode images 
%   meanFlnc - vector of surface fluence used to normalize images
%          NOTE: can be left empty if no normalization is to be performed. 
%   lambdas - vector of wavelengths (NOTE: used for burst titles only)
%   nPlt - specific bursts numbers to plot 
% 
% OUTPUT:
%   hp - figure handle of PA image.
% 
% Author: Vinoin Devpaul Vincely (2024)
% Edit History:
%   11/03/2024 - Included a detailed helper function. 
%   01/20/2024 - Transducer information for depth/lateral dimension
%   10/01/2020 - File Created! 


% Get depth information 
if size(BmodeImg,2) == 325
    [depth] = get_VSX_depth([15 size(BmodeImg,1)], 4, 15.625e6); 
    latDim = [-0.6 0.6]; climit = [6 8]; 
    climPA = [5.5 8];%-0.5;
elseif size(BmodeImg,2) == 129
    % L7-4 transducer
    [depth] = get_VSX_depth([0 220], 2, 5.208e6); 
    latDim = [-18 18]./10; climit = [7 9]; 
    climPA = [6 9];%-1;
else 
    error('This transducer has not been defined in "visualize_VSX_sPAI_data"!');
end 


figure(1); 
imagesc(log10(mean(BmodeImg,3))); colorbar; colormap(gray); clim([7 9]); 
% xlim([60 100]); ylim([20 100]); 
    tk.Yvals = linspace(1,size(BmodeImg,1),10); 
    tk.Ydpth = round(linspace(depth(1), depth(end), length(tk.Yvals)), 2);
    yticklabels(tk.Ydpth); yticks(tk.Yvals);
    tk.Xvals = linspace(1,size(BmodeImg,2),10); 
    tk.Xdpth = round(linspace(latDim(1), latDim(end), length(tk.Xvals)), 2);
    xticklabels(tk.Xdpth); xticks(tk.Xvals); 
    colormap(gray); colorbar; clim(climit); 
xlabel('Later Direction [cm]'); ylabel('Depth [cm]');
% hb = axes; 

figure(4); 
for c = 1:length(nPlt)
    hp(c) = subplot(1,length(nPlt),c); 
    if isempty(meanFlnc)
        imagesc(log10(mean(PA{nPlt(c)},3)));
    else 
        imagesc(log10(mean(PA{nPlt(c)},3)./meanFlnc(c))); % Fluence normalized PA image 
    end 
    colorbar; title(lambdas(nPlt(c)));     
    clim(climPA); % CLIMIT PA
%     xlim([60 100]); ylim([20 100]); 
        tk.Yvals = linspace(1,size(BmodeImg,1),10); 
        tk.Ydpth = round(linspace(depth(1), depth(end), length(tk.Yvals)), 2);
        yticklabels(tk.Ydpth); yticks(tk.Yvals);
        tk.Xvals = linspace(1,size(BmodeImg,2),10); 
        tk.Xdpth = round(linspace(latDim(1), latDim(end), length(tk.Xvals)), 2);
        xticklabels(tk.Xdpth); xticks(tk.Xvals);
    xlabel('Later Direction [cm]'); ylabel('Depth [cm]');
end 
colormap(hot); 
