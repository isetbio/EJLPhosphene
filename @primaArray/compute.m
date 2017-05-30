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



numberElectrodesX = floor(primaArray.width/primaArray.pixelWidth)+4;
numberElectrodesY = floor(primaArray.width/primaArray.pixelWidth)+4;
numberElectrodes = numberElectrodesX*numberElectrodesY;

activationWindow = ceil(size(movieInput,1)/numberElectrodesX);
%% Compute electrode activations from image

% Get the full image/movie from the identity outersegment
fullStimulus = movieInput;

% Find electrode activations by taking mean within window
for xPos = 1:numberElectrodesX
    for yPos = 1:numberElectrodesY
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
        electrodeStimulusL = squeeze(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:floor(imageCoordX1+activationWindow/2),:,:));
        electrodeStimulusR = squeeze(fullStimulus(imageCoordY1:imageCoordY2,floor(imageCoordX1+activationWindow/2)+1:imageCoordX2,:,:));
        % primaArray.activation(xPos,yPos,frame) = min([ mean(electrodeStimulus(:,1:floor(sizeES(2)/2))) mean(electrodeStimulus(:,ceil(sizeES(2)/2):sizeES(2)))]);
        primaArray.activation(xPos,yPos,:) = min([mean(RGB2XWFormat(electrodeStimulusL)); mean(RGB2XWFormat(electrodeStimulusL))]);
        
        else
            
            
        electrodeStimulusL = squeeze(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:floor(imageCoordX2),:,:));
        electrodeStimulusR = squeeze(fullStimulus(imageCoordY1:imageCoordY2,floor(imageCoordX1)+1:imageCoordX2,:,:));
        % primaArray.activation(xPos,yPos,frame) = min([ mean(electrodeStimulus(:,1:floor(sizeES(2)/2))) mean(electrodeStimulus(:,ceil(sizeES(2)/2):sizeES(2)))]);
        primaArray.activation(xPos,yPos,:) = min([mean(RGB2XWFormat(electrodeStimulusL)); mean(RGB2XWFormat(electrodeStimulusL))]);
        end
    end
end

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
clear bpMosaic


% Hack cone mosaic
cMosaicNS.pattern = 3*ones(size(primaArray.activationDS,1),size(primaArray.activationDS,2));
cMosaicNS.current = primaArray.activationDS;

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
        otherwise
            cMosaicNS.current = 1*primaArray.activationDSoff;
    end
    
    bpMosaic{cellTypeInd}.compute(cMosaicNS);
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
scaleFactor = 24;
innerRetina = scaleRF(innerRetina, scaleFactor);

innerRetina.compute(bpMosaic);

primaArray.innerRetina = innerRetina;
