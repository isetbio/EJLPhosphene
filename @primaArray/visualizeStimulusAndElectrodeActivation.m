function visualizeStimulusAndElectrodeActivation(primaArray, filename, electrodeCoords, fullStimulus, linearActivation, activation, activationDS, activationDSoff)

    % normalize to [0 .. 1]
    maxActivation = max([max(activation(:)) max(linearActivation(:))]);
    minActivation = min([min(activation(:)) min(linearActivation(:))]);
    maxStimulus = max(fullStimulus(:));
    minStimulus = min(fullStimulus(:));
    linearActivation = (linearActivation - minActivation)/(maxActivation-minActivation);
    activation = (activation - minActivation)/(maxActivation-minActivation);
    activationDS = (activationDS - minActivation)/(maxActivation-minActivation);
    activationDSoff = (activationDSoff - minActivation)/(maxActivation-minActivation);
    fullStimulus = (fullStimulus - minStimulus)/(maxStimulus-minStimulus);
    
    for k = 1:size(electrodeCoords,1)
        xAxis(k) = electrodeCoords(k,1).x;
    end
    for k = 1:size(electrodeCoords,2)
        yAxis(k) = electrodeCoords(1,k).y;
    end
    electrodeOutline.x = 0.5*(xAxis(2)-xAxis(1))*[-1 -1 1 1 -1];
    electrodeOutline.y = 0.5*(yAxis(2)-yAxis(1))*[-1 1 1 -1 -1];

    % Video writer
    videoOBJ = VideoWriter(sprintf('%s.mp4', filename), 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 5; 
    videoOBJ.Quality = 100;
    videoOBJ.open();
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 3, ...
       'heightMargin',  0.07, ...
       'widthMargin',    0.04, ...
       'leftMargin',     0.02, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.04);

    hFig = figure(1); clf; set(hFig, 'Color', [1 1 1], 'Position', [10 10 1200 800]);
    for f = 1:size(activation,3)
        subplot('Position', subplotPosVectors(1,1).v);
        imagesc(squeeze(fullStimulus(:,:,f)));
        hold on;
        for xPos = 1:size(electrodeCoords,1)
        for yPos = 1:size(electrodeCoords,2)
            plot(electrodeCoords(xPos,yPos).x, electrodeCoords(xPos,yPos).y, 'x', 'MarkerEdgeColor', electrodeCoords(xPos, yPos).rgb);
            plot(electrodeCoords(xPos,yPos).x + electrodeOutline.x, electrodeCoords(xPos,yPos).y + electrodeOutline.y, '-', 'Color', electrodeCoords(xPos, yPos).rgb)
        end
        end
        set(gca, 'CLim', [0 1], 'FontSize', 12);
        axis 'image'
        hold off
        title(sprintf('stimulus (frame %d) and electrode array', f), 'FontSize', 14);

        
        subplot('Position', subplotPosVectors(1,2).v);
        activationFrame = squeeze(linearActivation(:,:,f));
        imagesc(xAxis, yAxis, activationFrame);
        hold on
        for xPos = 1:size(electrodeCoords,1)
        for yPos = 1:size(electrodeCoords,2)
            plot(electrodeCoords(xPos,yPos).x + electrodeOutline.x, electrodeCoords(xPos,yPos).y + electrodeOutline.y, '-', 'Color', electrodeCoords(xPos, yPos).rgb)
        end
        end
        hold off;
        set(gca, 'CLim', [0 1], 'FontSize', 12);
        axis 'image'
        title('linear activation', 'FontSize', 14);
        
        
        subplot('Position', subplotPosVectors(1,3).v);
        activationFrame = squeeze(activation(:,:,f));
        imagesc(activationFrame);
        imagesc(xAxis, yAxis, activationFrame);
        hold on
        for xPos = 1:size(electrodeCoords,1)
        for yPos = 1:size(electrodeCoords,2)
            plot(electrodeCoords(xPos,yPos).x + electrodeOutline.x, electrodeCoords(xPos,yPos).y + electrodeOutline.y, '-', 'Color', electrodeCoords(xPos, yPos).rgb)
        end
        end
        hold off;
        set(gca, 'CLim', [0 1], 'FontSize', 12);
        axis 'image'
        title('min[e1,e2] activation', 'FontSize', 14);

        
        pos = subplotPosVectors(2,1).v;
        pos(1) = pos(1) + 0.15;
        subplot('Position', pos);
        activationFrame = squeeze(activationDS(:,:,f));
        imagesc(xAxis, yAxis, activationFrame);
        hold on
        for xPos = 1:size(electrodeCoords,1)
        for yPos = 1:size(electrodeCoords,2)
            plot(electrodeCoords(xPos,yPos).x + electrodeOutline.x, electrodeCoords(xPos,yPos).y + electrodeOutline.y, '-', 'Color', electrodeCoords(xPos, yPos).rgb)
        end
        end
        hold off;
        set(gca, 'CLim', [0 1], 'FontSize', 12);
        axis 'image'
        title('activationDS', 'FontSize', 14);
        
        pos = subplotPosVectors(2,2).v;
        pos(1) = pos(1) + 0.2;
        subplot('Position', pos);
        activationFrame = squeeze(activationDSoff(:,:,f));
        imagesc(xAxis, yAxis, activationFrame);
        hold on
        for xPos = 1:size(electrodeCoords,1)
        for yPos = 1:size(electrodeCoords,2)
            plot(electrodeCoords(xPos,yPos).x + electrodeOutline.x, electrodeCoords(xPos,yPos).y + electrodeOutline.y, '-', 'Color', electrodeCoords(xPos, yPos).rgb)
        end
        end
        hold off;
        set(gca, 'CLim', [0 1], 'FontSize', 12);
        axis 'image'
        title('activationDSoff', 'FontSize', 14);
        
        colormap(gray)
        drawnow
        if (f == 32)
            NicePlot.exportFigToPDF(sprintf('%sFrame%d.pdf', filename, f), hFig, 300)
        end
        videoOBJ.writeVideo(getframe(hFig));
    end
    videoOBJ.close();
end
