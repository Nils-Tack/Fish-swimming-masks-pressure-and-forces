function outlineFixed = repairOutline(outline,I,scalex,scaley)
% Select repair location and start repair on outline

%% Options
% Buffers for zoomed outline plots
% buff_lg = 50;
buff_sm = 20;

%% Show outline
figure
plot(outline(:,1),outline(:,2))
axis(bufferAxlims(outline(:,1),outline(:,2),10))
axis('equal')

figSize(1.2,1.2,gcf)
title('Click near section to repair.')

% Get axes around section to repair

[xax, yax] = ginput(1);
axlims_sm = [xax-buff_sm xax+buff_sm yax-buff_sm yax+buff_sm];

%% Replace points
fprintf('Replacing points.')
[outlineFixed,~,~] = replaceOutlinePoints(outline,I,scalex,scaley,axlims_sm);
