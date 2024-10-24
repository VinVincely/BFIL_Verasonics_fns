function [ROI_mask] = get_invivo_ROI_from_VSX(BmodeImg, PAImg)
% 
% function [ROI_mask] = get_invivo_ROI_from_VSX(Bmode, PA)
% 
% The objective of this funtion is to run a custom segmentation of the
% in vivo data collected with the verasonics system.

% Boilerplate checking 
if iscell(BmodeImg)
    error('Please input individual images!');
elseif size(PAImg,3) ~= 1 
    error('Please input individual images!');
end 

% Get depth information 
if size(BmodeImg,2) == 325
    [depth] = get_VSX_depth([15 size(BmodeImg,1)], 4, 15.625e6); 
    latDim = [-0.6 0.6]; climit = [6 8]; 
elseif size(BmodeImg,2) == 129
    [depth] = get_VSX_depth([0 220], 2, 5.208e6); 
    latDim = [-18 18]; climit = [7 9]; climit = [5.5 8.5];
else 
    % error('This transducer has not been defined in "get_invivo_ROI_from_VSX"!');
    [depth] = get_VSX_depth([0 220], 2, 5.208e6); 
    latDim = [-18 18]; climit = [7 9]; climit = [5.5 8.5];
end 

% plot data 
figure(1); 
imagesc(log10(BmodeImg)); xlabel('Lateral Distance [cm]'); ylabel('Depth [cm]');  
tk.Yvals = linspace(1,size(BmodeImg,1),10); 
tk.Ydpth = round(linspace(depth(1), depth(end), length(tk.Yvals)), 2);
yticklabels(tk.Ydpth); yticks(tk.Yvals);
tk.Xvals = linspace(1,size(BmodeImg,2),10); 
tk.Xdpth = round(linspace(latDim(1), latDim(end), length(tk.Xvals)), 2);
xticklabels(tk.Xdpth); xticks(tk.Xvals); 
colormap(gray); colorbar; clim(climit); 

figure(2); 
imagesc(log10(PAImg)); xlabel('Lateral Distance [cm]'); ylabel('Depth [cm]'); colormap(hot); 
tk.Yvals = linspace(1,size(BmodeImg,1),10); 
tk.Ydpth = round(linspace(depth(1), depth(end), length(tk.Yvals)), 2);
yticklabels(tk.Ydpth); yticks(tk.Yvals);
tk.Xvals = linspace(1,size(BmodeImg,2),10); 
tk.Xdpth = round(linspace(latDim(1), latDim(end), length(tk.Xvals)), 2);
xticklabels(tk.Xdpth); xticks(tk.Xvals);
% clim([7.5 9.6]);

% Draw ROI 
fprintf('Draw ROI using freehand...\n NOTE: Do not let-go of cursor until completion of drawing ROI!\n'); 
figure(1); ROI_info = drawfreehand; 
ROI_mask = createMask(ROI_info);





% ROI_OPEN = true; 
% x = 1;
% while ROI_OPEN
%     figure(1); line = drawline; 
%     line_pos = round(line.Position); 
%     lineX(x,:) = [line_pos(1,1) line_pos(2,1)]; 
%     lineY(x,:) = [line_pos(1,2) line_pos(2,2)]; 
%     if x == 1
%         start_pos = line_pos(1,:); 
%     end
%     if line_pos(2,1) == start_pos(1) && line_pos(2,2) == start_pos(2) 
%         ROI_OPEN = false;
%     end 
%     x = x + 1; 
% end 
% keyboard
%     if line_pos{1}
