function BWoutline = BWmasksFromOutlines2(imageSize,outline,scaleFactor)
BWoutline = flipud(poly2mask(round(outline(:,1)*scaleFactor),round(outline(:,2)*scaleFactor),imageSize(1),imageSize(2))); % rows,columns
end