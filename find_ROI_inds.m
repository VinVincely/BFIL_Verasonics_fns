function [objROI, bkROI] = find_ROI_inds(PA, Bmode, lims) 
% 
% function [objROI, bkROI] = find_ROI_inds(PA, Bmode, lims) 
% 
% The objective of this function is to find the indicies of rectangular
% "regions of interest" (ROI) for our signal and noise. The user will be
% prompted to draw each of these regions in sequence. The indicies will be
% saved in structures that are retured. 
% 
% INPUT:
%   PA    - PA image (Type: 2D array int/double)
%   Bmode - Bmode image (Type: 2D array int/double)
%   lims  - x, y and c limits for visualization (Type: struct, optional)
%       "x" - x axis limits (Type: 2-value vector int)
%       "y" - y axis limits (Type: 2-value vector int)
%       "c" - contrast axis limits (Type: 2-value vector int)
% 
% OUTPUT:
%   objROI - Indicies of signal ROI (Type: struct)
%   bkROI  - Indicies of background ROI (Type: struct)
%
% Author: Vinoin Devpaul Vincely (2024)

% Edit Log: 
%   07/03/2024 - Included a detailed helper funciton   

% ===== Input checks 
if isa(PA, 'double') 
    if size(PA,2) ~= 2
        error('"PA" must be an image of type double!'); 
    end    
else
    error('Please ensure "PA" is an image (2D double vector)!'); 
end 

if isa(Bmode, 'double') 
    if size(Bmode,2) ~= 2
        error('"Bmode" must be an image with two dimensions!');
    end    
else
    error('"Bmode" must be an image of type double!'); 
end 

if ~isa(lims, 'struct')
    error('"lim" must be a structure!'); 
end 

% ===== Limits included 




% keyboard

figure(1); imagesc(Bmode); xlabel('Channel no.'); ylabel('Depth [cm]');
colormap(gray);
% Fixing y axis 
depthAcq = size(Bmode,1)/2;
depth = get_VSX_depth([0 depthAcq], upsample, 5.208e6); 
ytick_vals = linspace(0,size(Bmode,1),10); 
depth = round(linspace(0, depth(end), length(ytick_vals)), 2);
yticklabels(depth); yticks(ytick_vals);

figure(2); imagesc(PA); xlabel('Channel no.'); ylabel('Depth [cm]');  
% Fixing y axis  
ytick_vals = linspace(0,size(PA,1),10); 
depth = round(linspace(0, depth(end), length(ytick_vals)), 2);
yticklabels(depth); yticks(ytick_vals); colormap(hot); 
%clim([0 5e7]); 
clim([5.5 7.5]);

if exist('lims', 'var')
    figure(1); ylim(lims); 
    figure(2); ylim(lims); colorbar;  %clim([0 1e8]); 
end

figure(2); 
fprintf('Please select the object ROI:\n'); 
r1 = drawrectangle('Label','obj','Color',[1 0 0]);
fprintf('Please select the background ROI:\n'); 
r2 = drawrectangle('Label','BK','Color',[1 0 0]);
    
objROI.chn = [ceil(r1.Vertices(1,1)):floor(r1.Vertices(3,1))]; 
objROI.dpth = [ceil(r1.Vertices(1,2)):floor(r1.Vertices(2,2))]; 
    
bkROI.chn = [ceil(r2.Vertices(1,1)):floor(r2.Vertices(3,1))]; 
bkROI.dpth = [ceil(r2.Vertices(1,2)):floor(r2.Vertices(2,2))]; 

% if isempty(objROI)
%     keyboard;
% end 