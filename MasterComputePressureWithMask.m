%% Compute pressure fields and forces from PIV
% Author: Nils Tack
% Compiled with Matlab 2020a
% Script used to prepare PIV data files, create masks/outlines from images,
% manually correct outlines (if necessary),compute pressure fields, thrust
% and drag forces acting on the body, and the Froude efficiency.

% Please refer to the following publication for details about pressure field 
% calculations from 2D velocity fields:
% Dabiri, J. O., Bose, S., Gemmell, B. J., Colin, S. P. and 
% Costello, J. H. (2014). An algorithm to estimate unsteady and quasi-steady 
% pressure fields from velocity field measurements. J. Exp. Biol. 217, 331–336.

% Before running the script, copy PIV data text files into the 'txt' folder
% found in 'Data'
% To make masks, also load corresponding image sequence into 'tif' folder
% found in 'Data'. Do not worry about frame increment yet. Copy all images based on
% start and end frames corresponding to the file number of first and last
% PIV text files.

% Required file format: image files (.jpg or .tif or .png), PIV data (.txt)

close all; clearvars

%% Main paths
path_main = 'D:\Documents\USF\Research\Research project\MATLAB tools\Matlab scripts\ComputePressureWithMask\Data'; % set path to 'Data' folder

% Other paths
path_PIVoriginal = fullfile(path_main,'PIVoriginal'); % raw PIV data
path_image = fullfile(path_main,'images'); % image sequence
path_PIVclean = fullfile(path_main,'PIVclean'); % velocity data without header
path_masksBW = fullfile(path_main,'masksBW'); % BW masks produced from original images
path_outlines = fullfile(path_main,'outlines'); % outlines produced from BW masks (.csv)
path_pressure = fullfile(path_main,'pressure'); % pressure data (.csv)
path_forces = fullfile(path_main,'forces'); % force data

% Set file type for image sequence (.jpg or .tif or .png)
imFileType = '.jpg';

%% OPTIONAL - calculate first and last frames of image sequence corresponding exactly to the first and last velocity field
calculateFirstLastFrame(400,430,2) % (first PIV file #,last PIV file #,increment set in DaVis to calculate velocity fields)
% To match exact frame # to PIV file number, this function perform the following calculation:(PIV file # - 1) * increment + 1

%% Select images based on the increment used for PIV data
incrPIV = 2;% same as increment originally set in DaVis to calculate velocity fields
selectImages(path_image,incrPIV,imFileType); % select images corresponding to each velocity field. Also renames files starting at 0001. File name becomes 'B".

%% OPTIONAL - Rotate images if necessary
% This step is completely optional. Though orientation will not affect the
% creation of masks and the calculation of the pressure fields, I recommend
% having the animal swimming up. This orientation is required for the
% computation of forces.
% To save time, it is best to orient all the original frames properly
% outside Matlab. The software used by Photron (PFV4) can be used to export
% rotated image sequences. WARNING! Rotating images using Microsoft’s
% rotate tool straight from the folder where all the images are found does
% NOT work.

opts.rotateSequence = 1;  % option to enable rotation of image sequence. 0 = no rotation, 1 = 90° rotation counter-clockwise

if opts.rotateSequence == 1
rotateImages(path_image,imFileType,180) % (folder containing images, file type, rotation angle); positive angle = CCW, negative angles = CW
end

%% Clean original PIV and export velocity data in the proper format (without column headers)
D_PIVoriginal = dir(fullfile(path_PIVoriginal,'*.txt')); % set directory; accepted file type .txt or .dat
filestub = 'B'; % Default file name
% 
fprintf('PIV files conversion...');
 for i=1:length(D_PIVoriginal)
      num = num2str(i,'%05d'); % create file number
      rf = strcat(filestub, num); % concatenate file number; filestub + file number
      A = importdata(fullfile(D_PIVoriginal(i).folder,D_PIVoriginal(i).name)); % read data
      csvwrite(fullfile(path_PIVclean,[rf,'.csv']),A.data(:,1:4)); % export clean data; accepted format .csv or .dat
      progressCount2(i,length(D_PIVoriginal)); % display export progress
  end
 fprintf('done\n');

%% Options for BW masks
opts.maskImage      = 0; % keep only area of interest (uses ginput)
opts.checkMasks     = 1; % check the quality of the masks as they get exported
opts.plotComparison = 1; % plots mask-making steps (mute automatically for mask export)
opts.inverseMask    = 0; % inverse the BW masks to make white areas black and vice-versa. Required if subject is dark and background white (brightfield PIV)

%% Test mask and image processing options
% Read 1 test image
i1 = 5; % test frame number
D_image = dir(fullfile(path_image,['*',imFileType])); %extract images from tif folder
temp = split(D_image(i1).name,'.');
filename = temp{1}; clear temp;
disp(filename)
I = importdata(quickfilepath(D_image(i1)));
I = I(:,:,1); % select only one of the three layers of the image

% Option to inverse the original grayscale image. makeMasks3 requires subject to be white on a black background -
% use this option if brightfield PIV.
if opts.inverseMask
    figure; hold on
    I = imcomplement(I);
    imshow(I)
    title('Inverse colors')
end

% Masking unnecessary parts of the image to isolate fish
if opts.maskImage
    [I,x_ginput,y_ginput] = maskMasked(I);
end

% Subplots to see if masking settings look good
% Open 'makeMask3' function to change parameters if necessary
% Parameters include masking unwanted areas like tank/flume edges or
% reflections, contrast and gamma corrections, smoothing options.
I4 = makeMask3(I, opts);

%% Export all the BW masks
opts.plotComparison = 0; % disable comparison figures

if opts.checkMasks % check individual masks as they get exported
    figTemp = figure;
end

fprintf('Exporting masks...');

for i=1:length(D_image)
    
    % Read original image data
    temp = split(D_image(i).name,'.');
    filename = temp{1}; clear temp;
    I = importdata(quickfilepath(D_image(i)));
    I = I(:,:,1); % Convert to grayscale
    
    % Option to inverse the original grayscale image. makeMasks3 requires subject to be white on a black background -
    % use this option if brightfield PIV.
    if opts.inverseMask
        I = imcomplement(I);
        imshow(I)
        title('Inverse colors')
    end

        % Masking unnecessary parts of the image to isolate fish
    if opts.maskImage
        ITemp = uint8(zeros(length(I(:,1)),length(I(1,:)))); % initiate new mask as black
        tempI = I(y_ginput(1):y_ginput(2),x_ginput(1):x_ginput(2));
        ITemp(y_ginput(1):y_ginput(2),x_ginput(1):x_ginput(2)) = tempI;
        I = ITemp;
    end
    
    % Make a binary mask
     I4 = makeMask3(I,opts);
        
    % Option to check masks during export
    if opts.checkMasks
       imshow(I4)
       title(sprintf('Frame %i',i))
       pause(0.05);
       clf
    end
        
    % Export mask
    filenameBW = sprintf('BW_%05g',i);
    imwrite(I4,fullfile(path_masksBW,[filenameBW,'.tif'])) % exports as tif     
end

close(figTemp)
fprintf('done\n');

%% Options for outlines
opts.smooth             = 3;    % 0 off, 1 moving average, 2 loess smoothing, 3 Savitzky-Golay smoothing
opts.smoothn            = 100;  % Smoothing parameter for opts.smooth
opts.visualizeExport    = 1;    % look at each outline during export to confirm the shape
pixelDilation           = 15;   % dilation factor (in px)

%% Set scale from PIV data to scale outlines
D_PIVclean = dir(fullfile(path_PIVclean,'*.csv')); % % directory for renamed PIV files; accepted file type .csv or .dat
D_BW = dir(fullfile(path_masksBW,'BW_*.tif')); % directory for BW masks

f = 1;  % select one image
BW = importdata(quickfilepath(D_BW(f))); % import mask
A = importdata(quickfilepath(D_PIVclean(f))); % import corresponding PIV file
X = A(:,1); % scale in mm
Y = A(:,2); % scale in mm

unqx=unique(X);
unqy=unique(Y);

% new scaling - matches DaVis output - vector placed at center of
% interrogation window. Accurate to +/- 1 px
scalex = linspace(min(unqx)-min(unqx),max(unqx)+min(unqx),size(BW,2))'; % increases width of x axis by 2*0.5 interrogation window
scaley = linspace(max(unqy)-max(unqy),min(unqy)+max(unqy),size(BW,1))'; % increases width of x axis by 2*0.5 interrogation window

% Axes limits and scale
axlims = [min(scalex), max(scalex), min(scaley), max(scaley)];
scale = size(BW,2)/(max(scalex)-min(scalex)); % px/mm

%% Test one outline
i1 = 1;
BW = flipud(importdata(quickfilepath(D_BW(i1)))); % import BW mask
I = importdata(quickfilepath(D_image(i1))); % import original image
I = I(:,:,1);

% Generate outline (unit is mm)
% WARNING - calculation of reliable pressure fields may require the dilation 
% of outlines to eliminate the boundary layer and to eliminar erroneous velocity 
% vectors potter in interrogation windows including the fish-water boundary. 
% Value (in px, depends greatly on image resolution)
outline = mask2outline(BW,scale,pixelDilation,opts); % function extracting the outline, (BW mask, scale in px/mm, dilate outline by n pixels, smoothing opts)

% Interpolate outline to generate evenly distributed points
N = 250; % number of points around outline
outline = curvspace(outline,N);

% plot original image with outline overlaid on it
figure; hold on
imagesc([axlims(1),axlims(2)],[axlims(3),axlims(4)],flipud(I)); colormap('gray')
plot(outline(:,1),outline(:,2),'.r','Linewidth',1.5);

% set(gca,'YDir','reverse')
axis equal

%% Export all the outlines
fprintf('Exporting outlines...');

if opts.visualizeExport
    figTemp = figure; 
end

for i=1:length(D_BW)
   
    % extract outline
    BW = flipud(importdata(quickfilepath(D_BW(i)))); %import mask
    outline = mask2outline(BW,scale,pixelDilation,opts); % function extracting the outline, (BW mask, scale in px/mm, dilate outline by n pixels, smoothing opts)
    outline = curvspace(outline,N);
    
    % Option to check the quality of the outlines
    if opts.visualizeExport
       I = importdata(quickfilepath(D_image(i))); % import original image
       I = I(:,:,1);
       
       hold on
       imagesc([min(scalex), max(scalex)],[min(scaley), max(scaley)],flipud(I)); colormap('gray')
       plot(outline(:,1),outline(:,2),'.r','Linewidth',1.5);
%        set(gca,'YDir','reverse')
       axis equal
       title(sprintf('Frame %i',i))
       pause(0.05);
       clf
    end
    
    % Export outline
    filename = sprintf('iface_%05g', i); % Number sequentially
    csvwrite(fullfile(path_outlines,[filename,'.csv']),outline)
    
end
close(figTemp)
fprintf('done\n');

%% OPTIONAL - Convert outlines to another unit (metric system)
opts.rescaleOutlines = 0; % option to rescale the outline to other metric unit. Set to 0 to prevent any accidental, unintended use of this function

if opts.rescaleOutlines == 1
d_mm = dir(fullfile(path_outlines,'*.csv'));
path_outline_mm = fullfile(path_outlines,'outlines_new_unit');

for i = 1:length(d_mm)
 temp = split(d_mm(i).name,'.'); % split file name
 filename = temp{1};

    B = []; % create empty temporary array
    A = readtable([d_mm(i).folder,filesep,d_mm(i).name]); % read x and y value for each outline file
    B(:,1) = A{:,1}./1000; % change value and operator to make any conversion in the metric system; i.e., /1000 converts from mm to m
    B(:,2) = A{:,2}./1000; % change value and operator to make any conversion in the metric system; i.e., /1000 converts from mm to m
    B = table(B(:,1),B(:,2));
    
poutline_mm_out = fullfile(path_outline_mm,[filename,'.csv']);% set export path
writetable(B,poutline_mm_out,'WriteVariableNames',0) % write table in the file format defined in the line above (.csv)
  
end
end

%% OPTIONAL - REPAIR OUTLINES MANUALLY IF NECESSARY
% Options
opts.exportOutlines = 1;   % Export outline files
opts.checkOutlines  = 1;   % Check existing outlines instead of making outline from mask
opts.keepTrackChange = 0;  % Keep track of which file has been checked

opts.repair         = 0;   % Repair mask before making outline
opts.smooth         = 2;   % 0 off, 1 moving average, 2 loess smoothing, 3 Savitzky-Golay smoothing
opts.smoothn        = 35;  % Smoothing parameter

opts.dilateMask     = 0;  % 0 = off, 1 = on. Option to dilate mask - not required if already performed upon creation of the masks
opts.dilatePx       = 0;  % number of pixels for dilation

%% Load files
fprintf('Loading files...');

% load BW masks
[masks_BW,filenums_BW] = getMasks(path_masksBW,2);
% mask_BW = masks_BW{1};

% load outlines
if opts.checkOutlines
    [outlines,filenums_out] = getMasks(path_outlines,0);
end
fprintf('done\n');
beep on; beep

%% Start repairing outlines
f1 = 1; % starting frame

% Loop through outlines
inter = 0;
for ff=f1:1:length(D_BW)
    temp = split(D_BW(ff).name,'.');
    filename = temp{1};
    filenumber = sscanf(temp{1},'BW_%d');
    disp(filename)
    
    BW = masks_BW{filenums_BW==filenumber};
    
    if opts.dilateMask
        SE = strel('disk',opts.dilatePx);
        BW = imdilate(BW,SE);
    end
    
    if opts.checkOutlines
        outline = outlines{filenums_out==filenumber};
    else
        outline = mask2outline(BW,scale,pixelDilation,opts);
        outline = curvspace(outline,N);
    end
    
    outlinePrev = outline;
    [I,~] = getImageForFrame(path_image,ff);
    
    inloop = 1; % Image loop
    changed = 0;
    while inloop
        close all
        showFrameWithOutline(I,outline,axlims)
        title([sprintf('%i - %s',ff,D_BW(ff).name),"Press 'a' to repair, 'k' to skip, 'z' to undo, 'return' to save & move to next, 's' to stop"],'interpreter','none')
       
        % Repair?
        keyChar = waitForKeyPress();
        switch keyChar
            case 'a' % Replace points
                outlinePrev = outline;
                outline = repairOutline(outline,I,scalex,scaley);
                changed = 1;
            case 'z' % Undo
                outline = outlinePrev; % Get previous outline
            case 's' % Stop
                inter = 1; % Interrupt outer loop
                break;
            case 'k' % Skip
                disp('Skipped')
                filename = [path_outlines,filesep,'Broken.csv'];
                dlmwrite(filename,filenumber,'-append')
                break; % Leave inner loop
            case 'p' % use outline from previous frame
                temp=split(D_BW(ff-1).name,',');
                filenumber_prev=sscanf(temp{1},'BW_%d');
                filename_prev=sprintf('%s%siface_%05i.csv', path_outlines,filesep,filenumber_prev);
                outline=importdata(filename_prev);
                changed=1;
            otherwise 
                if opts.keepTrackChange
                    filename = [path_outlines,filesep,'Checked.csv'];
                    dlmwrite(filename,filenumber,'-append')
                end
                break; % Leave inner loop
        end
    end
    
    if inter
        close all
        disp('Stopping')
        break
    end
    
% Save outlines with evenly spaced points
N=250; % number of points; change if needed
outline = curvspace(outline,N);
 
% Export outline
    if opts.exportOutlines
        if ~opts.checkOutlines || changed % If checking outlines, only export if outline has changed
            filename = sprintf('iface_%05g', filenumber); % Use original number
            fprintf('Exporting: %s\n',filename)
            csvwrite([path_outlines,filesep,filename,'.csv'],outline)
        end
    end
end
clf
fprintf('All masks repaired\n');

%% OPTIONAL - export corrected BW masks using the repaired outlines
opts.reExportBW = 0;
if opts.reExportBW == 1
    
D_outlines = dir(fullfile(path_outlines,'*.csv'));
imageSize = size(BW);

fprintf('Exporting BW masks from outlines...');
for i = 1:length(D_outlines)
    outline = importdata(quickfilepath(D_outlines(i))); % import outline
    BWoutline = BWmasksFromOutlines2(imageSize,outline,scale);
    
    % Export mask
    filenameBW = sprintf('BW_%05g',i);
    imwrite(BWoutline,fullfile(path_masksBW,[filenameBW,'.tif'])) % exports as tif
    
    progressCount2(i,length(D_outlines)); % display export progress
end
fprintf('done\n');
end

%% CALCULATE PRESSURE FIELDS
% Open parameters used by queen2 to run computation
% Change the required fields in the struct element
open('parameters2')

%% Run queen2
% Updated colormap for live pressure fields plot to enhance readability
% (see custom colormap lines 421 & 1090 for details)
queen2

%% Move pressure files from the PIVclean folder to pressure folder
D_preliminaryPress = dir([path_PIVclean,filesep,'press*.csv']);
fprintf('Moving press files...');

for id = 1:length(D_preliminaryPress)
    filepath_in = fullfile(path_PIVclean,D_preliminaryPress(id).name);
    filepath_out = fullfile(path_pressure,D_preliminaryPress(id).name);
    movefile(filepath_in, filepath_out);
end

fprintf('done\n');

%% COMPUTE FORCES ALONG THE BODY
% Options
opts.npoints         = 150;         % Number of points along the body to trace force vectors
opts.forceScale      = 3;          % Scale factor for plotting force vectors (for plotting purposes only, does not affect computation of forces)
opts.increment       = 1;           % Force calculation increment (from one frame to the next); generally stays at 1
opts.pressureDeltaT  = 0.002;       % Time interval between pressure fields in second (refer to steps above to determine the time between each pressure field)
opts.useRawImage     = 1;           % Option to plot raw image or mask (from outline); 0 = uses mask only (from outline), 1 = uses the raw image and overlays force vectors
opts.exportFigures   = 1;           % Option to export figures in the 'forces' folder
opts.swimmingSpeed   = 0.169401;    % Average swimming speed in meter per second of the animal assuming no acceleration

% Set directories
D_outlines = dir(fullfile(path_outlines,'*.csv')); % set directory (if not already set)
D_press = dir(fullfile(path_pressure,'*.csv')); % set directory for pressure files
D_image = dir(fullfile(path_image,['*',imFileType])); % set directory (if not already set)

% Test one frame
f1 = 126;
% options for plot: 'axial', 'lateral', 'both', 'all'
forcesAlongBody(D_outlines,D_press,D_image,opts,'both','test',f1);
% [summaryDataForce] = forcesAlongBody(D_outlines,D_press,D_image,opts,'axial','test',5); % for testing, add 'test' followed by frame of interest

%% Export all force figures and data
[summaryDataForce,allLocalForceData] = forcesAlongBody(D_outlines,D_press,D_image,opts,'axial');

%% Plot net force over time
timex = linspace(0,size(summaryDataForce,1)*opts.pressureDeltaT-opts.pressureDeltaT,size(summaryDataForce,1))';
figure('units','centimeters','Position',[1 1 15 15]); hold on
plot(timex,summaryDataForce(:,6),'-k','linewidth',2);
yline(0,'--k');
xlabel('Time (s)')
ylabel('Net thrust (N m^{-1})')

%% Calculate average Froude efficiency for 1 tail beat
startBeat = 12; % frame # where tail beat starts
endBeat = 119;   % frame # where tail beat ends
FroudePerBeat = mean(summaryDataForce(startBeat:endBeat,13))
