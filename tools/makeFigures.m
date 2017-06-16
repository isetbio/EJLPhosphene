function makeFigures(primaRecon, destinationFolder)


%%
cellTypeInd=2;
szBP = size(primaRecon.bpMosaic{cellTypeInd}.responseCenter);
micronsPerBP = 2;
figure; imagesc(micronsPerBP*([1:szBP(1)]-szBP(1)/2),micronsPerBP*([1:szBP(2)]-szBP(2)/2),(primaRecon.bpMosaic{cellTypeInd}.responseCenter(:,:,10)));
view(-90,90);
% 
% bpActivations = permute(primaRecon.bpMosaic{cellTypeInd}.responseCenter(:,:,10),[2 1]);
% figure; 
% [xi,yi] = meshgrid([1:szBP(1)],[1:szBP(2)]);
% % micronsPerBP*([1:szBP(1)]-szBP(1)/2),micronsPerBP*([1:szBP(2)]-szBP(2)/2)
% maxx = max(xi(:)); maxy = max(yi(:));
% scatter(micronsPerBP*(xi(:)-maxx/2),micronsPerBP*(yi(:)-maxy/2),2,bpActivations(:));
cellType = primaRecon.bpMosaic{cellTypeInd}.cellType;
title(sprintf('%s Bipolar Activations', cellType));
colormap gray
axis equal
axis off
% xlabel('Distance (m)'); ylabel('Distance (m)');
% set(gca,'fontsize',14);
print(gcf,'-dpng','/Users/james/Documents/MATLAB/EJLPhosphene/local/figures/june14/bipolar_offdiff.png')
%%

%%

% figure; imagesc(primaRecon.activation(:,:,10))
% 
% figure; imagesc(primaRecon.innerRetina.mosaic{1}.responseLinear(:,:,10))
%%
figure;
th = (0:1/6:1)'*2*pi;
scaleFactor = 1e6;
xh = scaleFactor*primaRecon.pixelWidth/2*cos(th);
yh = scaleFactor*primaRecon.pixelWidth/2*sin(th);

% % Plot electrode array
eaSize = size(primaRecon.center);
% figure;
hold on;
for j = 1:eaSize(1)
    for i = 1:eaSize(2)
        %         scatter(primaRecon.center(i,j,1),primaRecon.center(i,j,2));
        patch(xh+scaleFactor*primaRecon.center(i,j,1),yh+scaleFactor*primaRecon.center(i,j,2),primaRecon.activation(i,j,10))
    end
end
colormap gray
axis equal
axis off
title('Electrode Array Activations');
% xlabel('Distance (m)'); ylabel('Distance (m)');
set(gca,'fontsize',14);
print(gcf,'-dpng','/Users/james/Documents/MATLAB/EJLPhosphene/local/figures/june14/electrodeArray.png')
%%

% primaRecon.innerRetina.mosaic{1}.plot('mosaic')
primaRecon.innerRetina.mosaic{3}.plot('mosaicFill'); colormap gray;
% title('On Midget RGC Spikes');
title('On Midget RGC Spikes');
% axis equal
axis off
view(270,-90)
% print(gcf,'-dpng','/Users/james/Documents/MATLAB/EJLPhosphene/local/figures/june14/rgc_onmidget.png')