function BW2 = makeMaskTest(I, threshold)
msk = I > threshold; % value to only keep pixels brighter than threshold
    props = regionprops(logical(msk), 'Area');
    allAreas = sort([props.Area]);
    msk = bwareaopen(msk, 25);
    msk = imclearborder(msk);
    bw2 = bwareafilt(msk,1);
    BW2 = imfill(bw2,'holes');
end