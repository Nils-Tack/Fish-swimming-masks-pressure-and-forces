function [summaryData,allLocalForceData] = forcesAlongBody(D_outlines,D_press,D_tif,options,plottingOptions,test,testFrame)
npoints = options.npoints;
forceScale = options.forceScale;
increment = options.increment;
pressureDeltaT = options.pressureDeltaT;
useRawImage = options.useRawImage;
exportFigures = options.exportFigures;
swimmingSpeed = options.swimmingSpeed;

% various options to change appearance of quiver plots
arrowThickness = 1;                 % thickness of quiver arrows 
arrowHead = 0.15;                   % arrow head proportion relative to length of quiver
pullThrustCol = [0.45 0.68 0.8];    % color code of pull thrust force vectors
pushThrustCol = [0.8 0.5 0.5];      % color code of push thrust force vectors
pullDragCol = [0 0.2 0.5];          % color code of pull drag force vectors
pushDradCol = [0.5 0 0.1];          % color code of push drag force vectors

% Set figure and force data export path
if nargin < 6
tempPath = D_press.folder;
splitPath = split(tempPath,filesep);
PathForces = fullfile(cell2mat(join(splitPath(1:end-1),filesep)),'forces');
PathLocalForces  = fullfile(cell2mat(join(splitPath(1:end-1),filesep)),'forces','local force-power');
mkdir(PathLocalForces);
end

if exportFigures == 1 && nargin < 6
PathForcesFigures = fullfile(PathForces,'figures'); 
mkdir(PathForcesFigures) % make new directory for figures
end

% Set time increment between force fields to compute power and efficiency
deltat = increment*pressureDeltaT; % time between frames for force calculation

% Set scale from pressure files
I = importdata(quickfilepath(D_tif(1))); % import mask
C = importdata(quickfilepath(D_press(1))); % import corresponding pressure file
unqxm = unique(C(:,1));% scale in m
unqym = unique(C(:,2));% scale in m

% Matches PIV file scale but in m
scalexm = linspace(min(unqxm)-min(unqxm),max(unqxm)+min(unqxm),size(I,2))'; % increases width of x axis by 2*0.5 interrogation window
scaleym = linspace(max(unqym)-max(unqym),min(unqym)+max(unqym),size(I,1))'; % increases width of x axis by 2*0.5 interrogation window


if nargin > 5
    first = testFrame;
    last = testFrame;
else
    first = 1;
    last = length(D_press);
end

if nargin < 6
   if strcmp('all',plottingOptions)
        figure('units','centimeters','Position',[1 1 30 10]);
   else
        figure('units','centimeters','Position',[1 1 20 20]);
   end
end

advanceFr = 1;                      % Initiate variable to advance  index of force data to be saved
summaryData = zeros(last,13);       % Initiate summary data table to be output
allLocalForceData = cell(last,1);   % Initiate cell array to store local force and power data

for fr = first:increment:last
% Read files
    % Image
    if useRawImage == 1
    I = importdata(quickfilepath(D_tif(fr))); % Import image for plotting purposes
    I = I(:,:,1); % select only one of the three layers of the image
    end
    
    % Mask file
    blank = importdata(quickfilepath(D_outlines(fr))); % import outlines blanking coordinates
    blank = curvspace(blank,npoints); % number of points to re-define outline is wanted
    blank = alignPoints(blank); % always align array with anteriormost point of animal as first element
    boundx = blank(:,1)/1000;  % convert x coordinate from mm to meters
    boundy = blank(:,2)/1000;  % convert y coordinate from mm to meters
   
    % Pressure file
    press = importdata(quickfilepath(D_press(fr))); % import pressure data
    [Xpress,Ypress] = meshgrid(unique(press(:,1)),unique(press(:,2)));
    PP = griddata(press(:,1),press(:,2),press(:,7),Xpress,Ypress);

    
% Test plot
if nargin > 5
figure('units','centimeters','Position',[1 1 30 15]);
subplot(1,2,1);
hold on
if useRawImage == 1
imagesc([min(scalexm), max(scalexm)],[min(scaleym), max(scaleym)],flipud(I)); colormap('gray')
end
plot(boundx,boundy,'.r','Linewidth',1.5);
axis equal
title(sprintf('Outline - Frame %i',fr))

ax(1) = subplot(1,2,2);
axisMagnitude = 3;
hold on
contourf(Xpress,Ypress,PP,75,'edgecolor','none') % plot pressure fields
caxis([-axisMagnitude axisMagnitude]) % pressure field min and max magnitude
colormap(ax(1),customcolormap_preset('orange-white-purple'))    
patch(boundx,boundy,'k') % overlay outlines
hcb = colorbar;
hcb.Label.String = 'relative pressure (Pa)';
axis equal
title(sprintf('Pressure field - Frame %i',fr))
end

% calculate pressure at body boundary
surfpress = griddata(Xpress,Ypress,PP,boundx,boundy);

% Calculate length of each surface segment (equivalent to surface area per unit depth) between surface points. 
surfdx = boundx - circshift(boundx,1);   % compute x-distance between surface points
surfdy = boundy - circshift(boundy,1);   % compute y-distance between surface points
darea = sqrt(surfdx.^2 + surfdy.^2);     % compute surface area (per unit depth) between surface points; also equivalent to segment length

% Compute normal unit vectors (assuming that the animal is swimming UP)
surfunitnormx = -surfdy./darea; % compute x-component of vector normal to surface
surfunitnormy = surfdx./darea; % compute y-component of vector normal to surface

% Identify surface points where pressure is pushing/pulling in the same/oposite direction of swimming (axial component y)
indpospull = find(surfunitnormy < 0 & surfpress < 0);   % find surface points where low pressure is pulling animal forward (assuming animal is swimming UP)
indpospush = find(surfunitnormy > 0 & surfpress > 0);   % find surface points where high pressure is pushing animal forward (assuming animal is swimming UP)
indnegpull = find(surfunitnormy > 0 & surfpress < 0);   % find surface points where low pressure is pulling animal backward (assuming animal is swimming UP)
indnegpush = find(surfunitnormy < 0 & surfpress > 0);   % find surface points where high pressure is pushing animal backward (assuming animal is swimming UP)

% Identify surface points where pressure is pushing/pulling laterally (lateral component x)
indleftpull = find(surfunitnormx < 0 & surfpress < 0);   % find surface points where low pressure is pulling animal forward (assuming animal is swimming UP)
indleftpush = find(surfunitnormx < 0 & surfpress > 0);   % find surface points where high pressure is pushing animal forward (assuming animal is swimming UP)
indrightpull = find(surfunitnormx > 0 & surfpress < 0);   % find surface points where low pressure is pulling animal backward (assuming animal is swimming UP)
indrightpush = find(surfunitnormx > 0 & surfpress > 0);   % find surface points where high pressure is pushing animal backward (assuming animal is swimming UP)


%% Plot vectors for pressure according to the 4 categories above
if strcmp('axial',plottingOptions)
    if nargin > 5
    figure('units','centimeters','Position',[1 1 20 20]);
    end
    
    hold on
        if useRawImage == 0 
           patch(boundx,boundy,'k')
        else
           imagesc([min(scalexm), max(scalexm)],[min(scaleym), max(scaleym)],flipud(I)); colormap('gray')
        end

    % Quiver plot of forces along body
    % WARNING: use of 'forceScale' for plotting purposes only
    quiver(boundx(indpospull,:),boundy(indpospull,:),surfpress(indpospull,:).*zeros(size(surfunitnormx(indpospull,:))),surfpress(indpospull,:).*surfunitnormy(indpospull,:).*darea(indpospull,:)*forceScale,0,'Color',pullThrustCol,'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [0 1 1]
    quiver(boundx(indpospush,:),boundy(indpospush,:),surfpress(indpospush,:).*zeros(size(surfunitnormx(indpospush,:))),surfpress(indpospush,:).*surfunitnormy(indpospush,:).*darea(indpospush,:)*forceScale,0,'Color',pushThrustCol,'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [1 0.5 1]
    quiver(boundx(indnegpull,:),boundy(indnegpull,:),surfpress(indnegpull,:).*zeros(size(surfunitnormx(indnegpull,:))),surfpress(indnegpull,:).*surfunitnormy(indnegpull,:).*darea(indnegpull,:)*forceScale,0,'Color',pullDragCol,'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [0 0 1]
    quiver(boundx(indnegpush,:),boundy(indnegpush,:),surfpress(indnegpush,:).*zeros(size(surfunitnormx(indnegpush,:))),surfpress(indnegpush,:).*surfunitnormy(indnegpush,:).*darea(indnegpush,:)*forceScale,0,'Color',pushDradCol,'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [1 0 0]
    
    % quiver(boundx(indpospull,:),boundy(indpospull,:),-surfpress(indpospull,:).*zeros(size(surfunitnormx(indpospull,:))),-surfpress(indpospull,:).*surfunitnormy(indpospull,:).*darea(indpospull,:)*forceScale,0,'Color','#5d52ef','LineWidth',1,'MaxHeadSize',0.15);   % original color was [0 1 1]
    % quiver(boundx(indpospush,:),boundy(indpospush,:),-surfpress(indpospush,:).*zeros(size(surfunitnormx(indpospush,:))),-surfpress(indpospush,:).*surfunitnormy(indpospush,:).*darea(indpospush,:)*forceScale,0,'Color','#d36a13','LineWidth',1,'MaxHeadSize',0.15);   % original color was [1 0.5 1]
    % quiver(boundx(indnegpull,:),boundy(indnegpull,:),-surfpress(indnegpull,:).*zeros(size(surfunitnormx(indnegpull,:))),-surfpress(indnegpull,:).*surfunitnormy(indnegpull,:).*darea(indnegpull,:)*forceScale,0,'Color','#281199','LineWidth',1,'MaxHeadSize',0.15);   % original color was [0 0 1]
    % quiver(boundx(indnegpush,:),boundy(indnegpush,:),-surfpress(indnegpush,:).*zeros(size(surfunitnormx(indnegpush,:))),-surfpress(indnegpush,:).*surfunitnormy(indnegpush,:).*darea(indnegpush,:)*forceScale,0,'Color','#932a0d','LineWidth',1,'MaxHeadSize',0.15);   % original color was [1 0 0]

        if useRawImage == 0 
           % Add scale bar
           plot([0.002; 0.002+0.01], [0.010; 0.010], '-k', 'LineWidth', 3)
           text(0.002,0.013,'1 cm','Color','k')

           % Add scale arrow for force vector
           quiver(0.002,0.002,0.001*forceScale,0,'MaxHeadSize',0.2,'LineWidth',2,'Color','k','MaxHeadSize',0.5) % 0.1mN cm-1 = 0.01 N m-1
           text(0.002,0.005,'0.1 mN cm^{-1}','Color','k')
        else
           % Add scale bar
           plot([0.002; 0.002+0.01], [0.010; 0.010], '-w', 'LineWidth', 3)
           text(0.002,0.013,'1 cm','Color','w')

           % Add scale arrow for force vector
           quiver(0.002,0.002,0.001*forceScale,0,'MaxHeadSize',0.2,'LineWidth',2,'Color','w','MaxHeadSize',0.5) % 0.1mN cm-1 = 0.01 N m-1
           text(0.002,0.005,'0.1 mN cm^{-1}','Color','w')
        end

    title(sprintf('Axial force - frame %i',fr))
    axis off
    axis equal
    
elseif strcmp('lateral',plottingOptions)
        if nargin > 5
        figure('units','centimeters','Position',[1 1 20 20]);
        end
        
        hold on
            if useRawImage == 0 
               patch(boundx,boundy,'k')
            else
               imagesc([min(scalexm), max(scalexm)],[min(scaleym), max(scaleym)],flipud(I)); colormap('gray')
            end

    quiver(boundx(indleftpull,:),boundy(indleftpull,:),surfpress(indleftpull,:).*surfunitnormx(indleftpull,:).*darea(indleftpull,:)*forceScale,surfpress(indleftpull,:).*zeros(size(surfunitnormx(indleftpull,:))),0,'Color','#562689','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [0 1 1]
    quiver(boundx(indleftpush,:),boundy(indleftpush,:),surfpress(indleftpush,:).*surfunitnormx(indleftpush,:).*darea(indleftpush,:)*forceScale,surfpress(indleftpush,:).*zeros(size(surfunitnormx(indleftpush,:))),0,'Color','#b35807','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx(indrightpull,:),boundy(indrightpull,:),surfpress(indrightpull,:).*surfunitnormx(indrightpull,:).*darea(indrightpull,:)*forceScale,surfpress(indrightpull,:).*zeros(size(surfunitnormx(indrightpull,:))),0,'Color','#562689','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx(indrightpush,:),boundy(indrightpush,:),surfpress(indrightpush,:).*surfunitnormx(indrightpush,:).*darea(indrightpush,:)*forceScale,surfpress(indrightpush,:).*zeros(size(surfunitnormx(indrightpush,:))),0,'Color','#b35807','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);

    if useRawImage == 0 
           % Add scale bar
           plot([0.002; 0.002+0.01], [0.010; 0.010], '-k', 'LineWidth', 3)
           text(0.002,0.013,'1 cm','Color','k')

           % Add scale arrow for force vector
           quiver(0.002,0.002,0.001*forceScale,0,'MaxHeadSize',0.2,'LineWidth',2,'Color','k','MaxHeadSize',0.5) % 0.1mN cm-1 = 0.01 N m-1
           text(0.002,0.005,'0.1 mN cm^{-1}','Color','k')
        else
           % Add scale bar
           plot([0.002; 0.002+0.01], [0.010; 0.010], '-w', 'LineWidth', 3)
           text(0.002,0.013,'1 cm','Color','w')

           % Add scale arrow for force vector
           quiver(0.002,0.002,0.001*forceScale,0,'MaxHeadSize',0.2,'LineWidth',2,'Color','w','MaxHeadSize',0.5) % 0.1mN cm-1 = 0.01 N m-1
           text(0.002,0.005,'0.1 mN cm^{-1}','Color','w')
    end

    title(sprintf('Lateral force - frame %i',fr))
    axis off
    axis equal
 
elseif strcmp('both',plottingOptions)
    if nargin > 5    
    figure('units','centimeters','Position',[1 1 20 20]);
    end
    
        hold on
            if useRawImage == 0 
               patch(boundx,boundy,'k')
            else
               imagesc([min(scalexm), max(scalexm)],[min(scaleym), max(scaleym)],flipud(I)); colormap('gray')
            end

    quiver(boundx(indpospull,:),boundy(indpospull,:),forceScale*surfpress(indpospull,:).*surfunitnormx(indpospull,:).*darea(indpospull,:),forceScale*surfpress(indpospull,:).*surfunitnormy(indpospull,:).*darea(indpospull,:),0,'Color',[0.45 0.68 0.8],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx(indpospush,:),boundy(indpospush,:),forceScale*surfpress(indpospush,:).*surfunitnormx(indpospush,:).*darea(indpospush,:),forceScale*surfpress(indpospush,:).*surfunitnormy(indpospush,:).*darea(indpospush,:),0,'Color',[0.8 0.5 0.5],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx(indnegpull,:),boundy(indnegpull,:),forceScale*surfpress(indnegpull,:).*surfunitnormx(indnegpull,:).*darea(indnegpull,:),forceScale*surfpress(indnegpull,:).*surfunitnormy(indnegpull,:).*darea(indnegpull,:),0,'Color',[0 0.2 0.5],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx(indnegpush,:),boundy(indnegpush,:),forceScale*surfpress(indnegpush,:).*surfunitnormx(indnegpush,:).*darea(indnegpush,:),forceScale*surfpress(indnegpush,:).*surfunitnormy(indnegpush,:).*darea(indnegpush,:),0,'Color',[0.5 0 0.1],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead)

        if useRawImage == 0 
           % Add scale bar
           plot([0.002; 0.002+0.01], [0.010; 0.010], '-k', 'LineWidth', 3)
           text(0.002,0.013,'1 cm','Color','k')

           % Add scale arrow for force vector
           quiver(0.002,0.002,0.001*forceScale,0,'MaxHeadSize',0.2,'LineWidth',2,'Color','k','MaxHeadSize',0.5) % 0.1mN cm-1 = 0.01 N m-1
           text(0.002,0.005,'0.1 mN cm^{-1}','Color','k')
        else
           % Add scale bar
           plot([0.002; 0.002+0.01], [0.010; 0.010], '-w', 'LineWidth', 3)
           text(0.002,0.013,'1 cm','Color','w')

           % Add scale arrow for force vector
           quiver(0.002,0.002,0.001*forceScale,0,'MaxHeadSize',0.2,'LineWidth',2,'Color','w','MaxHeadSize',0.5) % 0.1mN cm-1 = 0.01 N m-1
           text(0.002,0.005,'0.1 mN cm^{-1}','Color','w')
        end

        title(sprintf('Force vectors - frame %i',fr))
        axis off
        axis equal
    
elseif strcmp('all',plottingOptions)
    if nargin > 5
    figure('units','centimeters','Position',[1 1 30 10]); 
    end
    
subplot(1,3,1)
    hold on
        if useRawImage == 0 
           patch(boundx,boundy,'k')
        else
           imagesc([min(scalexm), max(scalexm)],[min(scaleym), max(scaleym)],flipud(I)); colormap('gray')
        end

    % Quiver plot of forces along body
    % WARNING: use of 'forceScale' for plotting purposes only
    quiver(boundx(indpospull,:),boundy(indpospull,:),surfpress(indpospull,:).*zeros(size(surfunitnormx(indpospull,:))),surfpress(indpospull,:).*surfunitnormy(indpospull,:).*darea(indpospull,:)*forceScale,0,'Color',[0.45 0.68 0.8],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [0 1 1]
    quiver(boundx(indpospush,:),boundy(indpospush,:),surfpress(indpospush,:).*zeros(size(surfunitnormx(indpospush,:))),surfpress(indpospush,:).*surfunitnormy(indpospush,:).*darea(indpospush,:)*forceScale,0,'Color',[0.8 0.5 0.5],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [1 0.5 1]
    quiver(boundx(indnegpull,:),boundy(indnegpull,:),surfpress(indnegpull,:).*zeros(size(surfunitnormx(indnegpull,:))),surfpress(indnegpull,:).*surfunitnormy(indnegpull,:).*darea(indnegpull,:)*forceScale,0,'Color',[0 0.2 0.5],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [0 0 1]
    quiver(boundx(indnegpush,:),boundy(indnegpush,:),surfpress(indnegpush,:).*zeros(size(surfunitnormx(indnegpush,:))),surfpress(indnegpush,:).*surfunitnormy(indnegpush,:).*darea(indnegpush,:)*forceScale,0,'Color',[0.5 0 0.1],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [1 0 0]
    
    % quiver(boundx(indpospull,:),boundy(indpospull,:),-surfpress(indpospull,:).*zeros(size(surfunitnormx(indpospull,:))),-surfpress(indpospull,:).*surfunitnormy(indpospull,:).*darea(indpospull,:)*forceScale,0,'Color','#5d52ef','LineWidth',1,'MaxHeadSize',0.15);   % original color was [0 1 1]
    % quiver(boundx(indpospush,:),boundy(indpospush,:),-surfpress(indpospush,:).*zeros(size(surfunitnormx(indpospush,:))),-surfpress(indpospush,:).*surfunitnormy(indpospush,:).*darea(indpospush,:)*forceScale,0,'Color','#d36a13','LineWidth',1,'MaxHeadSize',0.15);   % original color was [1 0.5 1]
    % quiver(boundx(indnegpull,:),boundy(indnegpull,:),-surfpress(indnegpull,:).*zeros(size(surfunitnormx(indnegpull,:))),-surfpress(indnegpull,:).*surfunitnormy(indnegpull,:).*darea(indnegpull,:)*forceScale,0,'Color','#281199','LineWidth',1,'MaxHeadSize',0.15);   % original color was [0 0 1]
    % quiver(boundx(indnegpush,:),boundy(indnegpush,:),-surfpress(indnegpush,:).*zeros(size(surfunitnormx(indnegpush,:))),-surfpress(indnegpush,:).*surfunitnormy(indnegpush,:).*darea(indnegpush,:)*forceScale,0,'Color','#932a0d','LineWidth',1,'MaxHeadSize',0.15);   % original color was [1 0 0]

        if useRawImage == 0 
           % Add scale bar
           plot([0.002; 0.002+0.01], [0.010; 0.010], '-k', 'LineWidth', 3)
           text(0.002,0.013,'1 cm','Color','k')

           % Add scale arrow for force vector
           quiver(0.002,0.002,0.001*forceScale,0,'MaxHeadSize',0.2,'LineWidth',2,'Color','k','MaxHeadSize',0.5) % 0.1mN cm-1 = 0.01 N m-1
           text(0.002,0.005,'0.1 mN cm^{-1}','Color','k')
        else
           % Add scale bar
           plot([0.002; 0.002+0.01], [0.010; 0.010], '-w', 'LineWidth', 3)
           text(0.002,0.013,'1 cm','Color','w')

           % Add scale arrow for force vector
           quiver(0.002,0.002,0.001*forceScale,0,'MaxHeadSize',0.2,'LineWidth',2,'Color','w','MaxHeadSize',0.5) % 0.1mN cm-1 = 0.01 N m-1
           text(0.002,0.005,'0.1 mN cm^{-1}','Color','w')
        end

    title(sprintf('Axial force - frame %i',fr))
    axis off
    axis equal

subplot(1,3,2)
    hold on
            if useRawImage == 0 
               patch(boundx,boundy,'k')
            else
               imagesc([min(scalexm), max(scalexm)],[min(scaleym), max(scaleym)],flipud(I)); colormap('gray')
            end

    quiver(boundx(indleftpull,:),boundy(indleftpull,:),surfpress(indleftpull,:).*surfunitnormx(indleftpull,:).*darea(indleftpull,:)*forceScale,surfpress(indleftpull,:).*zeros(size(surfunitnormx(indleftpull,:))),0,'Color','#562689','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);   % original color was [0 1 1]
    quiver(boundx(indleftpush,:),boundy(indleftpush,:),surfpress(indleftpush,:).*surfunitnormx(indleftpush,:).*darea(indleftpush,:)*forceScale,surfpress(indleftpush,:).*zeros(size(surfunitnormx(indleftpush,:))),0,'Color','#b35807','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx(indrightpull,:),boundy(indrightpull,:),surfpress(indrightpull,:).*surfunitnormx(indrightpull,:).*darea(indrightpull,:)*forceScale,surfpress(indrightpull,:).*zeros(size(surfunitnormx(indrightpull,:))),0,'Color','#562689','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx(indrightpush,:),boundy(indrightpush,:),surfpress(indrightpush,:).*surfunitnormx(indrightpush,:).*darea(indrightpush,:)*forceScale,surfpress(indrightpush,:).*zeros(size(surfunitnormx(indrightpush,:))),0,'Color','#b35807','LineWidth',arrowThickness,'MaxHeadSize',arrowHead);

    if useRawImage == 0 
           % Add scale bar
           plot([0.002; 0.002+0.01], [0.010; 0.010], '-k', 'LineWidth', 3)
           text(0.002,0.013,'1 cm','Color','k')

           % Add scale arrow for force vector
           quiver(0.002,0.002,0.001*forceScale,0,'MaxHeadSize',0.2,'LineWidth',2,'Color','k','MaxHeadSize',0.5) % 0.1mN cm-1 = 0.01 N m-1
           text(0.002,0.005,'0.1 mN cm^{-1}','Color','k')
        else
           % Add scale bar
           plot([0.002; 0.002+0.01], [0.010; 0.010], '-w', 'LineWidth', 3)
           text(0.002,0.013,'1 cm','Color','w')

           % Add scale arrow for force vector
           quiver(0.002,0.002,0.001*forceScale,0,'MaxHeadSize',0.2,'LineWidth',2,'Color','w','MaxHeadSize',0.5) % 0.1mN cm-1 = 0.01 N m-1
           text(0.002,0.005,'0.1 mN cm^{-1}','Color','w')
    end

    title(sprintf('Lateral force - frame %i',fr))
    axis off
    axis equal

subplot(1,3,3)
    hold on
        if useRawImage == 0 
            patch(boundx,boundy,'k')
        else
            imagesc([min(scalexm), max(scalexm)],[min(scaleym), max(scaleym)],flipud(I)); colormap('gray')
        end

    quiver(boundx(indpospull,:),boundy(indpospull,:),forceScale*surfpress(indpospull,:).*surfunitnormx(indpospull,:).*darea(indpospull,:),forceScale*surfpress(indpospull,:).*surfunitnormy(indpospull,:).*darea(indpospull,:),0,'Color',[0.45 0.68 0.8],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx(indpospush,:),boundy(indpospush,:),forceScale*surfpress(indpospush,:).*surfunitnormx(indpospush,:).*darea(indpospush,:),forceScale*surfpress(indpospush,:).*surfunitnormy(indpospush,:).*darea(indpospush,:),0,'Color',[0.8 0.5 0.5],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx(indnegpull,:),boundy(indnegpull,:),forceScale*surfpress(indnegpull,:).*surfunitnormx(indnegpull,:).*darea(indnegpull,:),forceScale*surfpress(indnegpull,:).*surfunitnormy(indnegpull,:).*darea(indnegpull,:),0,'Color',[0 0.2 0.5],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead);
    quiver(boundx(indnegpush,:),boundy(indnegpush,:),forceScale*surfpress(indnegpush,:).*surfunitnormx(indnegpush,:).*darea(indnegpush,:),forceScale*surfpress(indnegpush,:).*surfunitnormy(indnegpush,:).*darea(indnegpush,:),0,'Color',[0.5 0 0.1],'LineWidth',arrowThickness,'MaxHeadSize',arrowHead)

        if useRawImage == 0 
           % Add scale bar
           plot([0.002; 0.002+0.01], [0.010; 0.010], '-k', 'LineWidth', 3)
           text(0.002,0.013,'1 cm','Color','k')

           % Add scale arrow for force vector
           quiver(0.002,0.002,0.001*forceScale,0,'MaxHeadSize',0.2,'LineWidth',2,'Color','k','MaxHeadSize',0.5) % 0.1mN cm-1 = 0.01 N m-1
           text(0.002,0.005,'0.1 mN cm^{-1}','Color','k')
        else
           % Add scale bar
           plot([0.002; 0.002+0.01], [0.010; 0.010], '-w', 'LineWidth', 3)
           text(0.002,0.013,'1 cm','Color','w')

           % Add scale arrow for force vector
           quiver(0.002,0.002,0.001*forceScale,0,'MaxHeadSize',0.2,'LineWidth',2,'Color','w','MaxHeadSize',0.5) % 0.1mN cm-1 = 0.01 N m-1
           text(0.002,0.005,'0.1 mN cm^{-1}','Color','w')
        end

        title(sprintf('Force vectors - frame %i',fr))
        axis off
        axis equal
end

if nargin < 6
    if exportFigures == 1
    %Export paths and options
%     FilepathSVG = fullfile(PathForcesFigures,sprintf('force%04d',fr)); % export .svg files (vector files that can be reworked in Illustrator)
    FilepathJPG = fullfile(PathForcesFigures,filesep,sprintf('force%04d',fr)); % export .jpg files
    % print('-depsc2', '-painters', FilepathSVG) %export eps
    %print('-dpng', Filepath) %export png
%     print('-dsvg', '-painters',FilepathSVG) %exports to svg (preferred)
    print(FilepathJPG,'-djpeg')
    %print('-dpdf', '-painters', '-bestfit', Filepath) %export pdf
    %print('-dtiff', filepath) %export tif
    end
end

if nargin < 6
   pause(0.1) % leaves time to visualize the plot. WARNING: increases the duration of the export
   hold off
   clf % clear figure to conserve memory
end
%% Calculate force acting on the body
% Axial force
totypospull = sum(darea(indpospull,:).*surfpress(indpospull,:).*surfunitnormy(indpospull,:));  % total forward pull force in this frame
totypospush = sum(darea(indpospush,:).*surfpress(indpospush,:).*surfunitnormy(indpospush,:));  % total forward push force in this frame
totynegpull = sum(darea(indnegpull,:).*surfpress(indnegpull,:).*surfunitnormy(indnegpull,:));  % total backward pull force in this frame
totynegpush = sum(darea(indnegpush,:).*surfpress(indnegpush,:).*surfunitnormy(indnegpush,:));  % total backward push force in this frame
netyforce = nansum(darea(:,:).*surfpress(:,:).*surfunitnormy(:,:)); % net force in this frame

% Lateral force
totxleftpull = sum(darea(indleftpull,:).*surfpress(indleftpull,:).*surfunitnormx(indleftpull,:));  % total left pull force in this frame
totxleftpush = sum(darea(indleftpush,:).*surfpress(indleftpush,:).*surfunitnormx(indleftpush,:));  % total left push force in this frame
totxrightpull = sum(darea(indrightpull,:).*surfpress(indrightpull,:).*surfunitnormx(indrightpull,:));  % total right pull force in this frame
totxrightpush = sum(darea(indrightpush,:).*surfpress(indrightpush,:).*surfunitnormx(indrightpush,:));  % total right push force in this frame
netxforce = nansum(darea(:,:).*surfpress(:,:).*surfunitnormx(:,:)); % net force in this frame       

%% Compute power and efficiency
if fr > 1 % cannot compute power of first frame because we need to know the position of the body in frame n-1
    % Import previous blanking outline file
    prevblank = importdata(quickfilepath(D_outlines(fr-increment))); % import previous outline blanking coordinates
    prevblank = curvspace(prevblank,npoints);
    prevblank = alignPoints(prevblank); % always align array with anteriormost point of animal as first element 
    prevboundx = prevblank(:,1)/1000;  % convert x coordinate from mm to meters
    prevboundy = prevblank(:,2)/1000;  % convert y coordinate from mm to meters
   
    %calculate local speed of body surface
    surfvelu = abs((boundx-prevboundx)./deltat); % x-direction
    surfvelv = abs((boundy-prevboundy)./deltat); % y-direction
    
    % calculate local force of animal on fluid
    localforcex = darea.*surfpress.*surfunitnormx; % x-direction
    localforcey = darea.*surfpress.*surfunitnormy; % y-direction
    
    % Calculate local power exerted by animal
    locallatpower = abs(surfvelu.*localforcex);   % in lateral direction (assuming animal swims in Y-DIRECTION)
    localaxipower = abs(surfvelv.*localforcey);   % in axial direction (assuming animal swims in Y-DIRECTION)
    
    % Calculate total power
    latpower = nansum(locallatpower,1); % lateral 
    axipower = nansum(localaxipower,1); % axial
    
    % Calculate total power exerted due to negative/positive suction/push pressure
    latpullpower = nansum(locallatpower(find(surfpress<0),1)); % total lateral power exerted due to low pressure suction
    latpushpower = nansum(locallatpower(find(surfpress>0),1)); % total lateral power exerted due to high pressure pushing
    axipullpower = nansum(localaxipower(find(surfpress<0),1)); % total axial power exerted due to low pressure suction
    axipushpower = nansum(localaxipower(find(surfpress>0),1)); % total axial power exerted due to high pressure pushing
    
    % Calculate fish swiming efficiency
    %swimmingSpeed = 0.060675; % average animal swimming speed in meter
    FroudeEfficiency = ((totypospull+totypospush)*swimmingSpeed)/(((totypospull+totypospush)*swimmingSpeed)+latpower);
    
else
    latpower = 0;
    axipower = 0;
    localforcex = darea.*surfpress.*surfunitnormx; % x-direction
    localforcey = darea.*surfpress.*surfunitnormy; % y-direction
    locallatpower = zeros(npoints,1);
    localaxipower = zeros(npoints,1);
    latpullpower = 0;
    latpushpower = 0;
    axipullpower = 0;
    axipushpower = 0;
    FroudeEfficiency = 0;
end

% Store summary force and power data for each frame in a new variable
% frame #, total pull thrust, total push thrust, total pull drag, total
% push drag, net force, total axial power, total lateral power, total axial
% pull power, total axial push power, total lateral pull power,  total lat
% push power, Froude efficiency
summaryData(advanceFr,:) = [fr,totypospull,totypospush,totynegpull,totynegpush,netyforce,axipower,latpower,axipullpower,axipushpower,latpullpower,latpushpower,FroudeEfficiency];

% Store local force and power along the body in cell array
% x coordinates (m), y coordinates (m), local lateral force component (N m-1), local axial
% force component (N m-1), local lateral power (W m-1), local axial power (W m-1)
localForcePowerData = [boundx,boundy,localforcex,localforcey,locallatpower,localaxipower]; % concatenate data
allLocalForceData(advanceFr) = {localForcePowerData};

% Export local force and power along the body (1 file per frame)
if nargin < 6
filename = sprintf('local_%04d.csv',fr); % file name and .csv file type
dlmwrite(fullfile(PathLocalForces,filename),localForcePowerData);
end

advanceFr = advanceFr+1; % advance index of summaryData by 1    
end

if nargin < 6
% Export final summary table
dataTable = array2table(summaryData);
dataTable.Properties.VariableNames = {'frame_#','total_pull_thrust_N.m-1','total_push_thrust_N.m-1','total_pull_drag_N.m-1','total_push_drag_N.m-1','net_force_N.m-1','total_axial_power_W.m-1','total_lateral_power_W.m-1','total_axial_pull_power_W.m-1','total_axial_push_power_W.m-1','total_lateral_pull_power_W.m-1','total_lateral_push_power_W.m-1','Froude efficiency'}; % set the name of each variable
% dataTable_out = fullfile(PathForces,'Forces and power summary.csv');
writetable(dataTable,fullfile(PathForces,'Forces and power summary.csv')) % write table in the file format defined in the line above (.csv)
fprintf('Summary table exported\n');
fprintf('Local force and power data exported\n');
close % close empty remaining figure
end

end