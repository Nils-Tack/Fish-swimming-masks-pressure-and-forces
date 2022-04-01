function progressCount2(i,last)
% lineLength = 0; % initialize number of characters to remove for each iteration
persistent lineLength;
if isempty(lineLength)
    lineLength = 0;
end

progress = round((i*100)/last); % calculate progress based on length of for loop
fprintf(repmat('\b',1,lineLength))% delete previous percentage in command line based on number of characters previously added
lineLength = fprintf('%i%%',progress); % displays progress and returns number of characters used
if i == last
    fprintf(repmat('\b',1,lineLength))
    lineLength = [];
end
end