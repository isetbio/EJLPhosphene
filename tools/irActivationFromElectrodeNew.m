function innerRetinaInput = irActivationFromElectrodeNew(innerRetina, electrodeArray, retinalPatchWidth, metersPerPixel, nTileRows, nTileCols, mosaicOffset, params, offFlag)

numberElectrodesX = floor(retinalPatchWidth/electrodeArray.width)+4;
numberElectrodesY = floor(retinalPatchWidth/electrodeArray.width)+4;
numberElectrodes = numberElectrodesX*numberElectrodesY;

innerRetinaInput = zeros([size(innerRetina.mosaic{1}.cellLocation),size(electrodeArray.activationDS,3)]);

% for frame = 1:params.nSteps
    for mosaicInd = 1:length(innerRetina.mosaic)
        
        if offFlag(mosaicInd)
            electrodeArrayCopy.activationDS = electrodeArray.activationDSoff;
        else
            electrodeArrayCopy.activationDS = electrodeArray.activationDS;
        end
        [xc, yc] = size(innerRetina.mosaic{mosaicInd}.cellLocation);
        for xind = 1:xc
            for yind = 1:yc
                % Find closest electrode
                % Electrode centers
                electrodeCenter = reshape(electrodeArray.center,[numberElectrodes,2]);
                % RGC center
                rgcCenterMosaic = innerRetina.mosaic{mosaicInd}.cellLocation{xind,yind}*metersPerPixel;
                [mosaicSubRow, mosaicSubCol] = ind2sub([nTileRows nTileCols],mosaicInd);
                rgcCenter = rgcCenterMosaic;% + squeeze(mosaicOffset(mosaicSubRow,mosaicSubCol,:))';
                % Find distance
                centerDistanceCoords = repmat(rgcCenter,[numberElectrodes,1]) - electrodeCenter;
                % Weight electrode activation by Gaussian according to distancej
                centerDistance = bsxfun(@(x,y)sqrt(x.^2+y.^2),centerDistanceCoords(:,1),centerDistanceCoords(:,2));
                centerDistanceRS = reshape(centerDistance,[numberElectrodesX,numberElectrodesY]);
                [centerDistanceOrd,centerDistanceInd] = sort(centerDistance,'ascend');
                
                [minDistance, minDistanceInd] = min(centerDistance);
                [xmin,ymin] = ind2sub([numberElectrodesX,numberElectrodesY],minDistanceInd);
%                 minXY(yind,xind,frame,mosaicInd,:) = [xmin ymin];
%                 innerRetinaInput(yind,xind,frame,mosaicInd) = electrodeArray.activation(xmin,ymin,frame)*exp(-minDistance/2e-4);
              
%                 innerRetinaInput(yind,xind,1:size(electrodeArrayCopy.activationDS(xmin,ymin,:),3),mosaicInd) = electrodeArrayCopy.activationDS(xmin,ymin,:)*exp(-minDistance/2e-4);

%                 %%
%                 lambdactr = 0;
%                 for lambda = 10.^([0:-.5:-8])
%                     
%                  innerRetinaInput = zeros([size(innerRetina.mosaic{1}.cellLocation),size(electrodeArray.activationDS,3)]);
                
                numberNearbyElectrodes = 7; % how many nearby electrodes over which to sum; 1 is closest, 7 is closest and hex pattern around
                for xind2 = 1:numberElectrodesX
                    for yind2 = 1:numberElectrodesY
                        if centerDistanceRS(xind2,yind2) <= centerDistanceOrd(numberNearbyElectrodes)
%                         innerRetinaInput(yind,xind,frame,mosaicInd) = ...
%                             innerRetinaInput(yind,xind,frame,mosaicInd) + ...
%                             electrodeArray.activation(xind2,yind2,frame)*exp(-centerDistanceRS(xind2,yind2)/2e-4);
                         innerRetinaInput(yind,xind,1:size(electrodeArrayCopy.activationDS(xmin,ymin,:),3),mosaicInd) = ...
                             innerRetinaInput(yind,xind,1:size(electrodeArrayCopy.activationDS(xmin,ymin,:),3),mosaicInd)+...
                             electrodeArrayCopy.activationDS(xind2,yind2,:)*exp(-centerDistanceRS(xind2,yind2)/.8e-5);
                
                        end
                    end
                end
%                 lambdactr =lambdactr+1;
%                outv(lambdactr) = innerRetinaInput(1,1,1,1)
%                 
%                 end
                %%
            end
        end
    end
% end
ph=1;
innerRetinaInput = innerRetinaInput./10;%(numberElectrodesX*numberElectrodesY);

% figure; imagesc(squeeze(minXY(:,:,10,1,1)))