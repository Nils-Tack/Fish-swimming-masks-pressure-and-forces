function rotateImages(path_tif,imFileType,angle)
D_tif = dir(fullfile(path_tif,['*',imFileType]));
fprintf('Image rotation...');
for i = 1:length(D_tif)
    I=imread(D_tif(i).name);
    rotated=imrotate(I,angle); % rotate frames counterclockwise
    temp = split(D_tif(i).name,imFileType);
    filename=fullfile(path_tif,[temp{1},imFileType]);
    imwrite(rotated,filename)
    progressCount2(i,length(D_tif)); % display export progress
end
fprintf('done\n');
end