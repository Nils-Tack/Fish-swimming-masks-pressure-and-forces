function [I3,x,y] = maskMasked(BW)
    figure;
    imshow(BW)
    
    I3 = uint8(zeros(length(BW(:,1)),length(BW(1,:)))); % initiate new mask as black
    
    [x,y] = ginput(2); % select two points to form the diagonal of the rectangle to be preserved
    x = fix(x); % prevents warning with next statements
    y = fix(y);% prevents warning with next statements

    tempI3 = BW(y(1):y(2),x(1):x(2));
    I3(y(1):y(2),x(1):x(2)) = tempI3;
    
end