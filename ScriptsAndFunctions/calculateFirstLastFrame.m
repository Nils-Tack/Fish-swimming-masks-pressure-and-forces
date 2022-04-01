function calculateFirstLastFrame(firstPIV,lastPIV,incr)
firstImage = (firstPIV-1)*incr+1;   % calculate frame number of image corresponding to first velocity field in the sequence
lastImage = (lastPIV-1)*incr+1;     % calculate frame number of image corresponding to last velocity field in the sequence

fprintf('Start of sequence: PIV field #%i --> image #%i\n',firstPIV,firstImage);
fprintf('End of sequence: PIV field #%i --> image #%i\n',lastPIV,lastImage);
end