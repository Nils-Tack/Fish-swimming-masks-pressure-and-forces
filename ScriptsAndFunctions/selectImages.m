function selectImages(path_tif,incr,imFileType)
% Copy tif files into a new folder based on set increment
fprintf('Image files extraction...');

TifFiles = dir([path_tif,filesep,'*',imFileType]); %path of folder for image files

filepath_Incr = [path_tif,filesep,'Incr=',num2str(incr)]; % create temporary folder
mkdir(filepath_Incr)

% copy images with increment into temporary folder
for i=1:incr:length(TifFiles)  
    filepath_SetIncr = [path_tif,filesep,TifFiles(i).name];
    copyfile(filepath_SetIncr,filepath_Incr);
end

% Delete the files not stored in the Incr=# folder and then moves the Incr=# files back into the main Tif folder
% Delete old files

for i=1:length(TifFiles)  % delete all the original files (no longer interesting)
    filepath_DeleteFile = [path_tif,filesep,TifFiles(i).name];
    delete(filepath_DeleteFile);
end

% Move new files back to Tif file
TifIncrFiles = dir([filepath_Incr,filesep,'*',imFileType]); %path of folder for Incr=# files (was jpg)
TifFilesBack = path_tif;
for i=1:length(TifIncrFiles)
    filepath_FileIncr = [filepath_Incr,filesep,TifIncrFiles(i).name];
    movefile(filepath_FileIncr,TifFilesBack);
end

% Delete empty Incr folder
rmdir(filepath_Incr);

% Re-number all the "tif" files
TifFiles = dir([path_tif,filesep,'*',imFileType]); % path of folder for image files

filestub = 'B';% rename file if necessary
for id = 1:length(TifFiles)
    num = num2str(id,'%05d');
    rf = strcat(filestub, num, imFileType);
    filepath_in = [path_tif,filesep,TifFiles(id).name];
    filepath_out = [path_tif,filesep,rf];
    movefile(filepath_in, filepath_out);
end

fprintf('done\n');

end