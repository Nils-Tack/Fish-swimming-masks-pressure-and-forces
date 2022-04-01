function renumberPIVclean(path_PIV)
PIVFiles = dir([path_PIV,filesep,'*.csv']);

fprintf('Renaming PIVclean...')
filestub = 'B-'; % be careful, matlab cannot rename a file with its own name. Must rename it to something different first, and then rename it to B
for id = 1:length(PIVFiles)
    num = num2str(id,'%05d');
    rf = strcat(filestub, num, '.csv');
    filepath_in = [path_PIV,filesep,PIVFiles(id).name];
    filepath_out = [path_PIV,filesep,rf];
    movefile(filepath_in, filepath_out);
end

PIVFiles = dir([path_PIV,filesep,'*.csv']);
filestub = 'B'; 
for id = 1:length(PIVFiles)
    num = num2str(id,'%05d');
    rf = strcat(filestub, num, '.csv');
    filepath_in = [path_PIV,filesep,PIVFiles(id).name];
    filepath_out = [path_PIV,filesep,rf];
    movefile(filepath_in, filepath_out);
end
fprintf('done\n');

end