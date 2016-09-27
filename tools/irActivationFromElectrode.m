function innerRetinaInput = irActivationFromElectrode(innerRetina, electrodeArray, retinalPatchWidth, metersPerPixel, nTileRows, nTileCols, mosaicOffset, params, offFlag)

if offFlag
    electrodeArrayCopy.activationDS = electrodeArray.activationDSoff;
else
    electrodeArrayCopy.activationDS = electrodeArray.activationDS;
end

numberElectrodesX = floor(retinalPatchWidth/electrodeArray.width)+4;
numberElectrodesY = floor(retinalPatchWidth/electrodeArray.width)+0;
numberElectrodes = numberElectrodesX*numberElectrodesY;

for frame = 1:params.nSteps
    for mosaicInd = 1:length(innerRetina.mosaic)
        [xc, yc] = size(innerRetina.mosaic{mosaicInd}.cellLocation);
        for xind = 1:xc
            for yind = 1:yc
                % Find closest electrode
                % Electrode centers
                electrodeCenter = reshape(electrodeArray.center,[numberElectrodes,2]);
                % RGC center
                rgcCenterMosaic = innerRetina.mosaic{mosaicInd}.cellLocation{xind,yind}*metersPerPixel;
                [mosaicSubRow, mosaicSubCol] = ind2sub([nTileRows nTileCols],mosaicInd);
                rgcCenter = rgcCenterMosaic + squeeze(mosaicOffset(mosaicSubRow,mosaicSubCol,:))';
                % Find distance
                centerDistanceCoords = repmat(rgcCenter,[numberElectrodes,1]) - electrodeCenter;
                % Weight electrode activation by Gaussian according to distancej
                centerDistance = bsxfun(@(x,y)sqrt(x.^2+y.^2),centerDistanceCoords(:,1),centerDistanceCoords(:,2));
                centerDistanceRS = reshape(centerDistance,[numberElectrodesX,numberElectrodesY]);
                
                [minDistance, minDistanceInd] = min(centerDistance);
                [xmin,ymin] = ind2sub([numberElectrodesX,numberElectrodesY],minDistanceInd);
                minXY(yind,xind,frame,mosaicInd,:) = [xmin ymin];
%                 innerRetinaInput(yind,xind,frame,mosaicInd) = electrodeArray.activation(xmin,ymin,frame)*exp(-minDistance/2e-4);
                innerRetinaInput(yind,xind,frame,mosaicInd) = electrodeArrayCopy.activationDS(xmin,ymin,frame)*exp(-minDistance/2e-4);
                
%                 for xind2 = 1:numberElectrodesX
%                     for yind2 = 1:numberElectrodesY
%                         innerRetinaInput(yind,xind,frame,mosaicInd) = ...
%                             innerRetinaInput(yind,xind,frame,mosaicInd) + ...
%                             electrodeArray.activation(xind2,yind2,frame)*exp(-centerDistanceRS(xind2,yind2)/2e-4);
%                     end
%                 end

            end
        end
    end
end

innerRetinaInput = innerRetinaInput./10;%(numberElectrodesX*numberElectrodesY);

% figure; imagesc(squeeze(minXY(:,:,10,1,1)))