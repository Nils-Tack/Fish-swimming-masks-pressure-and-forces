function showFrameWithOutline(I,outline,axlims)
[xmin, xmax, ymin, ymax] = deal(axlims(1),axlims(2),axlims(3),axlims(4));
figure; hold on
imagesc([xmin,xmax],[ymin,ymax],flipud(I)); colormap('gray')
set(gca,'YDir','normal')
plot(outline(:,1),outline(:,2),'.')
figSize(1.2,1.2,gcf);
axis([xmin xmax ymin ymax])
axis('square')
