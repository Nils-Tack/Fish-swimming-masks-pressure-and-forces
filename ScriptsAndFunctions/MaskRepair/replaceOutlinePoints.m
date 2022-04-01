function [outlineRepaired,newsect,nn] = replaceOutlinePoints(outline,I,scalex,scaley,axlims)
% Used to remove stirring rod from beginning of accel. sequences.
% Needed or does addMaskPts work?
%% Test
% disp('Testing')
% outlinebak = outline;
% axlims = axlims_sm;

%% Choose points
figure; hold on
imagesc(scalex,scaley,flipud(I)); colormap('gray')
plot(outline(:,1),outline(:,2),'r.')
% [~,dir,~] = rotationDir(mask); % Check direction
set(gca,'YDir','normal')
axis(axlims)
%axis('equal')
figSize(1.2,1.2,gcf)
% title('Select stray points. Close window when finished.')
title("Repair desired section of outline. Hit 'return' to finish.")
[xd,yd] = ginputc('ShowPoints',true,'Color','r'); % Use red crosshairs

%% Find nearest neighbors
nn = knnsearch(outline,[xd yd]);
% if nn(1,1) < nn(end,1)
%     nn = flipud(nn);  
% end
%% Is this the tail?
if max(nn) - min(nn) > size(outline(:,1))/2
    disp('tail')
    isTail = 1;
    % Shift by a quarter turn
    outline = circshift(outline,round(size(outline(:,1))/4));
    % Redo nearest neighbor search
    nn = knnsearch(outline,[xd yd]);
else
    isTail = 0;
end

%% Get new section
newsect = [xd,yd];
% Replace ends with nns
newsect(1,:) = outline(nn(1),:);
newsect(end,:) = outline(nn(end),:);

%% Check new section
% figure; hold on
% plot(newsect(1,1),newsect(1,2),'ro')
% plot(newsect(:,1),newsect(:,2))


%% Interpolate
newsect = interp2path(unique(newsect,'stable','row'),7,'pchip',1);

%% Replace points

if nn(1) < nn(end)
% Split mask into two sections, excluding new section
mask_start = outline(1:nn(1),:);
mask_end = outline(nn(end):end,:);

% Fill in section
outlineRepaired = [mask_start; newsect(2:end-1,:); mask_end]; %original 
    
else
mask_start = outline(1:nn(end),:);% was 1
mask_end = outline(nn(1):end,:);% was end

% Fill in section
outlineRepaired = [mask_start; flipud(newsect(2:end-1,:)); mask_end]; %original 

end

outlineRepaired = unique(outlineRepaired,'stable','row');

% Close path
% outlineRepaired(end+1,:) = outlineRepaired(1,:);
% outlineRepaired(end,:) = [];

if isTail
    % Shift back
    outlineRepaired = circshift(outlineRepaired,round(size(outline(:,1))/4));
end



% %% Replace points
% % Split mask into two sections, excluding new section
% mask_start = outline(1:nn(1),:);
% mask_end = outline(nn(end):end,:);
% 
% 
% % Fill in section
% outlineRepaired = [mask_start; newsect(2:end-1,:); mask_end]; %original 
%     
% 
% 
% % Close path
% % outlineRepaired(end+1,:) = outlineRepaired(1,:);
% outlineRepaired(end,:) = [];
% 
% if isTail
%     % Shift back
%     outlineRepaired = circshift(outlineRepaired,round(size(outline(:,1))/4));
% end

%% Show new mask
% figure; hold on
% imagesc(scalex,scaley,flipud(I)); colormap('gray') % show image
% plot(outlineRepaired(:,1),outlineRepaired(:,2),'r.-')
% % plot(maskNoSect(:,1),maskNoSect(:,2),'r.')
% plot(newsect(:,1),newsect(:,2),'bo');
% axis(axlims)
% figSize(1.5,1.5,gcf)

