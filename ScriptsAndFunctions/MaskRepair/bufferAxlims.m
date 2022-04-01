function axlims = bufferAxlims(xdata,ydata,buffer)
% Get axis limits fit to data with buffer on all sides
    axlims = [min(xdata)-buffer, max(xdata)+buffer, min(ydata)-buffer, max(ydata)+buffer];
end