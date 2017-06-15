function compute(primaArray, movieInput, varargin)
%COMPUTE - a function of primaArray for computing the object.


%%

p = inputParser;
addRequired(p, 'movieInput');
p.KeepUnmatched = true;

p.parse(movieInput,varargin{:});
movieInput = p.Results.movieInput;

%% Cone Mosaic
% Generate the cone mosaic for retinal properties
% The absorptions and photocurrent will not be used

coneParams.fov = primaArray.fov;

iStimNS = ieStimulusMovieCMosaic(rand(100,100,1),coneParams);
cMosaicNS = iStimNS.cMosaic;

%% Prima



numberElectrodesX = floor(primaArray.width/primaArray.pixelWidth)+1;
numberElectrodesY = floor(primaArray.width/primaArray.pixelWidth)+1;
numberElectrodes = numberElectrodesX*numberElectrodesY;

activationWindow = ceil(size(movieInput,1)/numberElectrodesX);
%% Compute electrode activations from image

% Get the full image/movie from the identity outersegment
fullStimulus = movieInput;

% Find electrode activations by taking mean within window
% imElec = zeros(100,100);

electrodeAtten = fspecial('Gaussian', round(activationWindow(1)), activationWindow(1)/2);
electrodeAtten = electrodeAtten./max(electrodeAtten(:));
electrodeAttenTemporal = repmat(electrodeAtten,[1 1 size(fullStimulus,3)]);
for xPos = 1:numberElectrodesX
    for yPos = 1:numberElectrodesY
%         [xPos yPos]
        % Xcoords of window for stimulus
        imageCoordX1 = (activationWindow)*(xPos-1)+1;
        imageCoordX2 = (activationWindow)*(xPos);
        
        % Ycoords of window for stimulus
        imageCoordY1 = (activationWindow)*(yPos-1)+1;
        imageCoordY2 = (activationWindow)*(yPos);
        
        
        if imageCoordX1 < 1; imageCoordX1 = 1; end;
        if imageCoordY1 < 1; imageCoordY1 = 1; end;
        
        if imageCoordX1 > size(fullStimulus,2); imageCoordX1 = size(fullStimulus,2); end;
        if imageCoordY1 > size(fullStimulus,1); imageCoordY1 = size(fullStimulus,1); end;
        if imageCoordX2 > size(fullStimulus,2); imageCoordX2 = size(fullStimulus,2); end;
        if imageCoordY2 > size(fullStimulus,1); imageCoordY2 = size(fullStimulus,1); end;
        % Pull out piece of stimulus and take mean
        electrodeStimulus = squeeze(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2,:,:));
        % primaArray.activation(xPos,yPos,:) = mean(RGB2XWFormat(electrodeStimulus));                
        
        % Implement the local electrode min([e1,e2]) nonlinearity
        sizeES = size(electrodeStimulus);
        if imageCoordX1 < (size(fullStimulus,2)-activationWindow/2)
%             electrodeAttenL = electrodeAttenTemporal(1:length(imageCoordY1:imageCoordY2), 1:length(imageCoordX1:floor(imageCoordX1+activationWindow/2)),:);
%             electrodeAttenR = electrodeAttenTemporal(1:length(imageCoordY1:imageCoordY2),end-length(floor(imageCoordX1+activationWindow/2)+1:imageCoordX2)+1:end,:);
            
            electrodeStimulusL = squeeze(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:floor(imageCoordX1+activationWindow/2),:,:));
            electrodeStimulusR = squeeze(fullStimulus(imageCoordY1:imageCoordY2,floor(imageCoordX1+activationWindow/2)+1:imageCoordX2,:,:));
            % primaArray.activation(xPos,yPos,frame) = min([ mean(electrodeStimulus(:,1:floor(sizeES(2)/2))) mean(electrodeStimulus(:,ceil(sizeES(2)/2):sizeES(2)))]);
            % primaArray.activation(xPos,yPos,:) = min([mean(RGB2XWFormat(electrodeAttenL.*electrodeStimulusL)); mean(RGB2XWFormat(electrodeAttenR.*electrodeStimulusR))]);
            primaArray.activation(xPos,yPos,:) = min([mean(RGB2XWFormat(1.*electrodeStimulusL)); mean(RGB2XWFormat(1.*electrodeStimulusR))]);
            
            %         imElec(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2) = (xPos*yPos)*ones(size(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2)));
        else
            
%             electrodeAttenL = electrodeAttenTemporal(1:length(imageCoordY1:imageCoordY2), 1:length(imageCoordX1:floor(imageCoordX2)),:);
%             electrodeAttenR = electrodeAttenTemporal(1:length(imageCoordY1:imageCoordY2),end-length(floor(imageCoordX1)+1:imageCoordX2)+1:end,:);
  
            electrodeStimulusL = squeeze(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:floor(imageCoordX2),:,:));
            electrodeStimulusRsq = squeeze(fullStimulus(imageCoordY1:imageCoordY2,floor(imageCoordX1)+1:imageCoordX2,:,:)); 
            szESR = size(((electrodeStimulusRsq)));
            electrodeStimulusR = zeros(szESR(1),1,szESR(2));
            electrodeStimulusR(1:szESR(1),:,1:szESR(2)) = electrodeStimulusRsq;
            % primaArray.activation(xPos,yPos,frame) = min([ mean(electrodeStimulus(:,1:floor(sizeES(2)/2))) mean(electrodeStimulus(:,ceil(sizeES(2)/2):sizeES(2)))]);
            % primaArray.activation(xPos,yPos,:) = min([mean(RGB2XWFormat(electrodeAttenL.*electrodeStimulusL)); mean(RGB2XWFormat(electrodeAttenR.*electrodeStimulusR))]);
            primaArray.activation(xPos,yPos,:) = min([mean(RGB2XWFormat(1.*electrodeStimulusL)); mean(RGB2XWFormat(1.*electrodeStimulusR))]);
            
            %         imElec(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2) = (xPos*yPos)*ones(size(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2)));
        end
    end
end
% figure; imagesc(imElec);
%% Add X Hz spiking of stimulus
% primaArray.pulseFreq = pulseFreq;
% Right now the electrode sampling is at 0.01 s = 100 Hz
% Downsample to get 5 Hz
szAct = size(primaArray.activation);
primaArray.activationDS = zeros(szAct);
for iSample = 1:size(fullStimulus,3)
    if mod(iSample,100/primaArray.pulseFreq)==0
        primaArray.activationDS(:,:,iSample) = primaArray.activation(:,:,iSample);
        primaArray.activationDSoff(:,:,iSample) = 1-primaArray.activation(:,:,iSample);
    end
end

eaRS = reshape(primaArray.activation,[szAct(1)*szAct(2),szAct(3)]);
eaDSRS = reshape(primaArray.activationDS,[szAct(1)*szAct(2),szAct(3)]);

%% Bipolar
% The general strategy is to take our electrode activations, which will be
% a relatively coarse grid, and set their activations into the current
% field of the cone mosaic. The bipolar responses to the electrode
% activations are then calculated at this low spatial resolution to save
% processing time. Then the electrode activations are scaled up to the full
% size of the bipolar mosaic.

clear bpMosaic


coneSize = size(cMosaicNS.pattern);
patchSizeMicrons = [cMosaicNS.height cMosaicNS.width];
% Hack cone mosaic
cMosaicNS.pattern = 3*ones(size(primaArray.activationDS,1),size(primaArray.activationDS,2));
% cMosaicNS.current = primaArray.activationDS;



cellType = {'ondiffuse','offdiffuse','onmidget','offmidget','onSBC'};
for cellTypeInd = 1:4
    clear bpParams
    bpParams.cellType = cellType{cellTypeInd};
    
    bpParams.ecc = primaArray.ecc;
    bpParams.rectifyType = 1;
    bpMosaic{cellTypeInd} = bipolar(cMosaicNS, bpParams);
    bpMosaic{cellTypeInd}.set('sRFcenter',1);
    bpMosaic{cellTypeInd}.set('sRFsurround',0);
    switch cellType{cellTypeInd}
        case {'ondiffuse','onparasol','onmidget'}            
            cMosaicNS.current = 1*primaArray.activationDS;
%             for fr = 1:size(primaArray.activationDS,3)
%                 cMosaicNS.current = [];
%                 cMosaicNS.current(:,:,fr) = imresize(primaArray.activationDS(:,:,fr),[coneSize],'method','nearest');
%             end
        otherwise
            cMosaicNS.current = 1*primaArray.activationDSoff;
%             for fr = 1:size(primaArray.activationDS,3)
%                 cMosaicNS.current(:,:,fr) = imresize(primaArray.activationDSoff(:,:,fr),[coneSize]);
%             end
    end
    
    bpMosaic{cellTypeInd}.compute(cMosaicNS);
    bpResponseCenterTemp = bpMosaic{cellTypeInd}.responseCenter;
    bpResponseSurroundTemp = bpMosaic{cellTypeInd}.responseSurround;
        
    
    cMosaicNS.pattern = 3*ones([coneSize]);
    cMosaicNS.current = [];
    cMosaicNS.current(:,:,1) = imresize(primaArray.activationDS(:,:,1),[coneSize],'method','nearest');
            
    bpMosaic{cellTypeInd} = bipolar(cMosaicNS, bpParams);
    
    bpMosaic{cellTypeInd}.set('sRFcenter',1);
    bpMosaic{cellTypeInd}.set('sRFsurround',0);
    
   
    bpSize = size(bpMosaic{cellTypeInd}.cellLocation);
    bipolarsPerMicron = bpSize(1)./(1e6*patchSizeMicrons(1));

    scaleFactor = 1e6 * bipolarsPerMicron *.98;
    electrodeStimCenter = zeros(coneSize(1),coneSize(2),size(primaArray.activationDS,3));
    electrodeStimSurround = zeros(coneSize(1),coneSize(2),size(primaArray.activationDS,3));
    centerOffset = [primaArray.height primaArray.width]/2;    
%     centerOffset = [2*primaArray.pixelWidth+primaArray.center(end,end,1)-primaArray.center(1,1,1) 2*primaArray.pixelWidth+primaArray.center(end,end,2)-primaArray.center(1,1,2)]/2;    
    
    activationWindow = ceil([coneSize(1),coneSize(2)]/numberElectrodesX);
    primaArray.spatialWeight = fspecial('Gaussian', round(1.25*activationWindow(1)), activationWindow(1)/2);
    primaArray.spatialWeight = primaArray.spatialWeight./max(primaArray.spatialWeight(:));
    szWeight = size(primaArray.spatialWeight);
    for ri = 1:size(primaArray.center,1)-0
        for ci = 1:size(primaArray.center,2)-0            
            
            centerCoords = scaleFactor*(centerOffset+[primaArray.center(ri,ci,1),primaArray.center(ri,ci,2)]);        
            rc = [ceil(-szWeight(1)/2+centerCoords(1)) : floor(szWeight(1)/2+centerCoords(1))]; 
            rcind1 = find(rc>0); rcind2 = find(rc<coneSize(1)); rcind =intersect(rcind1,rcind2); rc = rc(rcind);
            cc = [ceil(-szWeight(2)/2+centerCoords(2)) : floor(szWeight(2)/2+centerCoords(2))]; 
            ccind1 = find(cc>0); ccind2 = find(cc<coneSize(1)); ccind =intersect(ccind1,ccind2); ccind =intersect(ccind1,ccind2); cc = cc(ccind);
           weightRS = primaArray.spatialWeight(rcind,ccind);
%             electrodeStimCenter(rc,cc,:) = XW2RGBFormat(weightRS(:)*squeeze(bpResponseCenterTemp(ri,ci,:))',length(rc), length(cc));            
%             electrodeStimSurround(rc,cc,:) = XW2RGBFormat(weightRS(:)*squeeze(bpResponseSurroundTemp(ri,ci,:))',length(rc), length(cc));      
            electrodeStimCenter(rc,cc,:) = electrodeStimCenter(rc,cc,:)+ XW2RGBFormat(weightRS(:)*squeeze(bpResponseCenterTemp(ri,ci,:))',length(rc), length(cc));            
            electrodeStimSurround(rc,cc,:) = electrodeStimSurround(rc,cc,:) + XW2RGBFormat(weightRS(:)*squeeze(bpResponseSurroundTemp(ri,ci,:))',length(rc), length(cc));         
            
        end
    end
%     figure; imagesc(electrodeStimMask(:,:,10))
    
    bpMosaic{cellTypeInd}.set('responseCenter',electrodeStimCenter);
    bpMosaic{cellTypeInd}.set('responseSurround',electrodeStimSurround);
end

primaArray.bpMosaic = bpMosaic;

%% RGC
clear params rgcParams
params.eyeRadius = primaArray.ecc;
params.eyeAngle = 90;
innerRetina=ir(bpMosaic,params);
cellType = {'on parasol','off parasol','on midget','off midget'};

rgcParams.centerNoise = 0;
rgcParams.model = 'LNP';
rgcParams.ellipseParams = [1 1 0];  % Principle, minor and theta

rgcParams.type = cellType{1};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{2};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{3};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{4};
innerRetina.mosaicCreate(rgcParams);
% 
% for mosaicInd = 1:length(innerRetina.mosaic)
%     for ri = 1:size((innerRetina.mosaic{mosaicInd}.cellLocation),1);
%         for ci = 1:size((innerRetina.mosaic{mosaicInd}.cellLocation),2);
%             sRFcenterNew{ri,ci} = 16*innerRetina.mosaic{mosaicInd}.sRFcenter{ri,ci};
%             sRFsurroundNew{ri,ci} = 16*innerRetina.mosaic{mosaicInd}.sRFsurround{ri,ci};
%         end
%     end
%     innerRetina.mosaic{mosaicInd}.set('sRFcenter', sRFcenterNew);
%     innerRetina.mosaic{mosaicInd}.set('sRFsurround', sRFsurroundNew);
% end
scaleFactor = 1;
innerRetina = scaleRF(innerRetina, scaleFactor);

innerRetina.compute(bpMosaic);

primaArray.innerRetina = innerRetina;
