function corrOutline = alignPoints(outline)
% 'outline' must be a n x 2 array with outline(:,1) = x coordinates and outline(:,2) = y coordinates
% because outlines are smoothed or reshaped (number of points changed)
% when imported to compute forces along the body, it is likely that edges
% of the matrix shift. Because the order of the the points matter when
% computing power (power calculation assuming change in space of points
% fixed on the surface of the animal), it is highly recommended to fix the
% edges of the matrix. For simplicity, the anteriormost point is selected
% as the beginning of the matrix.
% WARNING: this works for animals whose anteriormost point is precisely 
% aligned with its anteroposterior axis (fish, gellyfish). 

boundx = outline(:,1);
boundy = outline(:,2);

k = find(boundy==max(boundy)); % find index of max value of y (uppermost point of the body)
newBoundx=circshift(boundx,-k(1)+1); % move array up; using k(1) just in case multiple points fall on the same y
newBoundy=circshift(boundy,-k(1)+1); % move array up

corrOutline = [newBoundx,newBoundy];

end