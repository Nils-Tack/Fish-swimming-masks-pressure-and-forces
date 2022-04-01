function outline = mask2outline(BW,scale,n,opts) 

% Dilate masks - useful to remove boundary layer for accurate pressure
% field calculation
BW = bwmorph(BW,'thicken',n);

% get the coordinates of boundaries of all the white areas
BWedge = bwboundaries(BW,'noholes');
BWoutline = BWedge{1};
outline = [BWoutline(:,2)/scale,BWoutline(:,1)/scale]; % set the scale for the outline

% prevents noise caused by smoothing that tends to keep the outline open
if opts.smooth
    [outline,ia,ic]=unique(outline(:,1:2),'stable','rows'); % remove duplicate rows
    smooth_error = 10; % arbitrary number of points taken from begining of raw outline to add to the end of outline to increase accuracy of location of last point output by smoothing; 10-15 is good
    outline = [outline;outline(1:smooth_error,:)];
    
% options to smooth the outline
switch opts.smooth
  
    case 1
        % Use moving average
        outline(:,1) = smooth(outline(:,1),opts.smoothn); % 3-5
        outline(:,2) = smooth(outline(:,2),opts.smoothn);
    case 2    
        % Use loess smoothing
        outline(:,1) = smooth(outline(:,1),opts.smoothn,'loess'); % 20
        outline(:,2) = smooth(outline(:,2),opts.smoothn,'loess');
    case 3
        % Use Savitzky-Golay smoothing
        outline(:,1) = smooth(outline(:,1),opts.smoothn,'sgolay'); % 20
        outline(:,2) = smooth(outline(:,2),opts.smoothn,'sgolay');
    case 0
        
end
 outline(end-smooth_error:end,:) = []; % removes the unnecessary points added before smoothing
end
end