function BWmasksFromOutlines(D_BW,D_outlines,path_masksBW,scale)
% create coordinate matrix of pixels in new mask
BW = flipud(importdata(quickfilepath(D_BW(1)))); % import 1 BW mask
BWx = 1:1:size(BW,1);
BWy = 1:1:size(BW,2);
[XX,YY] = meshgrid(unique(BWx),unique(BWy));

fprintf('Exporting BW masks from outlines...');
for i = 1:length(D_outlines)
    % Convert outline from mm to px
    outlineTest = importdata(quickfilepath(D_outlines(i)))*scale;
    
    % Find pixels in and on ouline (returns logical)
    in = flipud(inpolygon(XX,YY,outlineTest(:,1),outlineTest(:,2)));
    
    % Export mask
    filenameBW = sprintf('BW_%05g',i);
    imwrite(in,fullfile(path_masksBW,[filenameBW,'.tif'])) % exports as tif  
end
fprintf('done\n');
end