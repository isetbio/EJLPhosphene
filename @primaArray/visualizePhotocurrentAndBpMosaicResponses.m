function visualizePhotocurrentAndBpMosaicResponses(primaArray, filename, ...
    weights, photocurrentResponse, bpResponseCenter, bpResponseCenterFull)
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', size(primaArray.center,1), ...
       'colsNum', size(primaArray.center,2), ...
       'heightMargin',  0.025, ...
       'widthMargin',    0.025, ...
       'leftMargin',     0.02, ...
       'rightMargin',    0.02, ...
       'bottomMargin',   0.02, ...
       'topMargin',      0.01);
   
    hFig = figure(10); clf; set(hFig, 'Color', [1 1 1], 'Position', [10 10 1050 1050]);
    for ri = 1:size(primaArray.center,2)
        for ci = 1:size(primaArray.center,1)
            weightRS = weights{ri,ci};
            subplot('Position', subplotPosVectors(ri,ci).v);
            imagesc(weightRS);
            axis 'image';
            set(gca, 'FontSize', 12);
            colormap(gray);
            drawnow
        end
    end
    NicePlot.exportFigToPDF(sprintf('%s_weights.pdf', filename), hFig, 300);
    
    
    % normalize to [0 .. 1]
    photocurrentRange = [min(photocurrentResponse(:)) max(photocurrentResponse(:))];
    bpCenterRange = max(abs(bpResponseCenter(:))) *[-1 1];
    photocurrentRange
    bpCenterRange
    
    photocurrentResponse = (photocurrentResponse - photocurrentRange(1))/(photocurrentRange(2)-photocurrentRange(1));
    bpResponseCenter = (bpResponseCenter-bpCenterRange(1))/(bpCenterRange(2)-bpCenterRange(1));
    bpResponseCenterFull = (bpResponseCenterFull -bpCenterRange(1))/(bpCenterRange(2)-bpCenterRange(1));
    
    for k = 1:size(primaArray.electrodeCoords,1)
        xAxis(k) = primaArray.electrodeCoords(k,1).x;
    end
    for k = 1:size(primaArray.electrodeCoords,2)
        yAxis(k) = primaArray.electrodeCoords(1,k).y;
    end
    
    electrodeOutline.x = 0.5*(xAxis(2)-xAxis(1))*[-1 -1 1 1 -1];
    electrodeOutline.y = 0.5*(yAxis(2)-yAxis(1))*[-1 1 1 -1 -1];
    
    xAxisFull = 1:size(bpResponseCenterFull,2);
    yAxisFull = 1:size(bpResponseCenterFull,1);
    dX = primaArray.electrodeCoordsFull(1,2).x - primaArray.electrodeCoordsFull(1,1).x;
    dY = primaArray.electrodeCoordsFull(2,1).y - primaArray.electrodeCoordsFull(1,1).y;
    electrodeOutlineFull.x = 0.5*(dX)*[-1 -1 1 1 -1];
    electrodeOutlineFull.y = 0.5*(dY)*[-1 1 1 -1 -1];
    
    
    % Video writer
    videoOBJ = VideoWriter(fullfile(phospheneRootPath(), sprintf('%s.mp4', filename)), 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 5; 
    videoOBJ.Quality = 100;
    videoOBJ.open();
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 1, ...
       'colsNum', 3, ...
       'heightMargin',  0.07, ...
       'widthMargin',    0.06, ...
       'leftMargin',     0.02, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.04);
   
    ticks = 0:0.125:1.0;
    for k = 1:numel(ticks)
        photocurrentTickLabels{k} = sprintf('%2.1f', photocurrentRange(1)+(k-1)/(numel(ticks)-1)*(photocurrentRange(2)-photocurrentRange(1)));
        bpCenterTicks{k} = sprintf('%2.1f', bpCenterRange(1)+(k-1)/(numel(ticks)-1)*(bpCenterRange(2)-bpCenterRange(1)));
    end
    
    hFig = figure(1); clf; set(hFig, 'Color', [1 1 1], 'Position', [10 10 1100 350]);
    
    for f = 1:size(photocurrentResponse,3)
        subplot('Position', subplotPosVectors(1,1).v);
        imagesc(xAxis, yAxis, squeeze(photocurrentResponse(:,:,f)));
        hold on;
        plotElectrodes(primaArray.electrodeCoords, electrodeOutline, true);
        hold off;
        set(gca, 'CLim', [0 1],  'FontSize', 12);
        colorbar('Ticks', ticks, 'TickLabels', photocurrentTickLabels);
        axis 'image'
        title(sprintf('photocurrent/ \nelectrode activation'),  'FontSize', 14);
        
        subplot('Position', subplotPosVectors(1,2).v);
        imagesc(xAxis, yAxis, squeeze(bpResponseCenter(:,:,f)));
        hold on;
        plotElectrodes(primaArray.electrodeCoords, electrodeOutline, true);
        hold off;
        set(gca, 'CLim', [0 1],  'FontSize', 12);
        colorbar('Ticks', ticks, 'TickLabels', bpCenterTicks);
        axis 'image'
        title(sprintf('bpResponseCenter\n(single electrode)'),  'FontSize', 14);
        
        subplot('Position', subplotPosVectors(1,3).v);
        imagesc(xAxisFull, yAxisFull, squeeze(bpResponseCenterFull(:,:,f)));
        set(gca, 'CLim', [0 1],  'FontSize', 12);
        colorbar('Ticks', ticks, 'TickLabels', bpCenterTicks);
        axis 'image'
        hold on
        plotElectrodes(primaArray.electrodeCoordsFull, electrodeOutlineFull, true);
        hold off;
        title(sprintf('bpResponseCenter\n(all electrodes)'),  'FontSize', 14);
        
        colormap(gray)
        drawnow;
        
        if (f == 32)
            NicePlot.exportFigToPDF(sprintf('%sFrame%d.pdf', filename, f), hFig, 300);
        end
        videoOBJ.writeVideo(getframe(hFig)); 
    end    
    videoOBJ.close();
end

function plotElectrodes(electrodeCoords, electrodeOutline, showElectrodeCenter)
    for xPos = 1:size(electrodeCoords,1)
        for yPos = 1:size(electrodeCoords,2)
            if (showElectrodeCenter)
                plot(electrodeCoords(xPos,yPos).x, electrodeCoords(xPos,yPos).y, 'g.', 'MarkerSize', 22, 'Color', electrodeCoords(xPos, yPos).rgb);
            end
            %plot(electrodeCoords(xPos,yPos).x + electrodeOutline.x, electrodeCoords(xPos,yPos).y + electrodeOutline.y, '-', 'Color', electrodeCoords(xPos, yPos).rgb)
        end
    end
end
