function [I,filename] = getImageForFrame(pathTIF,f)
% Get image for frame f

D = dir([pathTIF,filesep,'*.jpg']);
temp = split(D(f).name,'.');
filename = temp{1}; clear temp;
disp(filename)
I = importdata(quickfilepath(D(f)));