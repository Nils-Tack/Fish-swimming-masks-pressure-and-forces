function filename = quickfilepath(dI)
% Take an item from a directory and return its fullpath
filename = fullfile(dI.folder,dI.name);