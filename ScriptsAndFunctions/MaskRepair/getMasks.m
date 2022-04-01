function [masks,filenums] = getMasks(pathMasks,BW)
% Read mask files
% BW: 0 - outlines, 1 - BW masks in csv format, 2 - BW masks in tif format, 3 - original tif files 

switch BW
    case 0
        D = dir([pathMasks,filesep,'iface_*.csv']);
    case 1
        D = dir([pathMasks,filesep,'BW_*.csv']);
    case 2
        D = dir([pathMasks,filesep,'BW_*.tif']);
    case 3
        D = dir([pathMasks,filesep,'*.tif']);
end

if isempty(D)
    disp('No files.')
end

masks = cell(length(D),1);
filenums = zeros(length(D),1);
for i=1:length(D)
    filename = quickfilepath(D(i));
%     disp(filename)
    
    switch BW
        case 0
            masks{i} = importdata(filename);
            filenums(i) = sscanf(D(i).name,['iface_', '%d', '.csv']);
        case 1
            masks{i} = importdata(filename);
            filenums(i) = sscanf(D(i).name,['BW_', '%d', '.csv']);
        case 2
            masks{i} = imread(filename);
            filenums(i) = sscanf(D(i).name,['BW_', '%d', '.tif']);
        case 3
            masks{i} = imread(filename);
            filenums(i) = sscanf(D(i).name,['B', '%d', '.tif']);
    end

end