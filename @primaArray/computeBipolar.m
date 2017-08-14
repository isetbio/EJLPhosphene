function primaArray = computeBipolar(primaArray, cMosaicNS)
%COMPUTEELECTRODE - a function of primaArray for computing the bipolar activations.
%
%   primaArray.computeBipolar(cMosaic);
%
% The general strategy is to take our electrode activations, which will be
% a relatively coarse grid, and set their activations into the current
% field of the dummy cone mosaic. The bipolar responses to the electrode
% activations are then calculated at this low spatial resolution to save
% processing time. Then the electrode activations are scaled up to the full
% size of the bipolar mosaic, and a spatial falloff in bipolar activation
% from the electrode is implemented.
% 
% 5/2017 JRG (c) isetbio team


%% Cone Mosaic
% Generate the dummy cone mosaic to get the properties of the retinal patch
% (size, etc.). The absorptions and photocurrent are not computed or used.

coneSize = size(cMosaicNS.pattern);
patchSizeMicrons = [cMosaicNS.height cMosaicNS.width];

cMosaicNS.pattern = 3*ones(size(primaArray.activationDS,1),size(primaArray.activationDS,2));
% cMosaicNS.current = primaArray.activationDS;

%% Compute bipolar response
% Since the bipolar mosaic must take a coneMosaic object as its input, we
% take the dummy cone mosaic and set the computed electrode activations
% into the 'current' field of the cone mosaic. This way, the bipolar
% current output is computed with the electrode activations as input.
% 
% The bipolar output is the result of a temporal convolution of the bipolar
% filter with the cone mosaic 'current' (which is actually the electrode
% activation signal). Since the size of the electrode array (e.g. 15x15xt) is
% much smaller than the bipolar array (256x256xt), it is fastest to compute the
% temporal convolution on the 15x15 electrode signal in the cone mosaic
% 'current', and then resize the output signal to the actual size of the
% bipolar mosaic.
% 
% There is also a spatial falloff implemented after the temporal
% convolution during the scaling to the full size of the bipolar array.

numberElectrodesX = size(primaArray.activation,1);
numberElectrodesY = size(primaArray.activation,2);

% Five different types of bipolar cells
cellType = {'ondiffuse','ondiffuse','onmidget','onmidget','onSBC'};


clear bpL

bpL = bipolarLayer(cMosaicNS);

% Make each type of bipolar mosaic
cellType = {'on diffuse','off diffuse','on midget','off midget'};%,'on sbc'};

% Stride isn't influencing yet.s
clear bpMosaicParams
bpMosaicParams.rectifyType = 1;  % Experiment with this
bpMosaicParams.spread  = 1;  % RF diameter w.r.t. input samples
bpMosaicParams.stride  = 1;  % RF diameter w.r.t. input samples
bpMosaicParams.spreadRatio  = 10;  % RF diameter w.r.t. input samples
bpMosaicParams.ampCenter = 1.5;%1.3;%1.5 _2
bpMosaicParams.ampSurround = .5;%1;%.5
% Maybe we need a bipolarLayer.compute that performs this loop
for cellTypeInd = 1:length(cellType)
    bpL.mosaic{cellTypeInd} = bipolarMosaic(cMosaicNS, cellType{cellTypeInd}, bpMosaicParams);
    % bpL.mosaic{cellTypeInd}.compute();
end


fprintf('In computeBipolar\n')

for cellTypeInd = 1%:4
    % Set the cone mosaic current field to either the positive or negative
    % activations depending on bipolar cell type
    switch cellType{cellTypeInd}
        case {'ondiffuse','onparasol','onmidget'}
            cMosaicNS.current = 1*primaArray.activationDS;
        otherwise
            cMosaicNS.current = 1*primaArray.activationDSoff;
    end
    currentTemp = cMosaicNS.current;
    
    % Apply the temporal filter by computing
    bpL.mosaic{cellTypeInd}.set('sRFcenter',1);
    bpL.mosaic{cellTypeInd}.set('sRFsurround',0);
    bpL.mosaic{cellTypeInd}.compute();
    bpResponseCenterTemp = bpL.mosaic{cellTypeInd}.responseCenter;
    bpResponseSurroundTemp = bpL.mosaic{cellTypeInd}.responseSurround;
    
    
    %%%%%%%%%%%%%
    % Resize the bipolar activation to the full mosaic and implement spatial falloff
    
    % Set all dummy cones to M cones
    cMosaicNS.pattern = 3*ones([coneSize]);
    
    % Set current field to upsampled version of the electrode activations
    % for the first time step here (this is done for the rest of the time
    % steps below).
    cMosaicNS.current = [];
    fprintf('imresize')
    size(primaArray.activationDS(:,:,1))
    coneSize
    
    cMosaicNS.current(:,:,1) = imresize(primaArray.activationDS(:,:,1),[coneSize],'method','nearest');
    
end
    
    
for cellTypeInd = 1:4
    % Generate the bpMosaic
    % clear bpL.mosaic{cellTypeInd}
    bpL.mosaic{cellTypeInd} = bipolarMosaic(cMosaicNS, cellType{cellTypeInd}, bpMosaicParams);
end
    
for cellTypeInd = 1%:4
    % Assume fovea and set relative size to cones as 1:1
%     bpL.mosaic{cellTypeInd}.set('sRFcenter',1);
%     bpL.mosaic{cellTypeInd}.set('sRFsurround',0);
    
    % Get the conversion factor from bipolars to microns
    bpSize = size(bpL.mosaic{cellTypeInd}.cellLocation);
    bipolarsPerMicron = bpSize(1)./(1e6*patchSizeMicrons(1));
    
    % Set the scale factor for the spatial attenuation
    scaleFactor = 1e6 * bipolarsPerMicron *.98;
    
    % Initialize full spatial scale response for bp
    bpResponseCenterFull = zeros(coneSize(1),coneSize(2),size(primaArray.activationDS,3));
    bpResponseSurroundFull = zeros(coneSize(1),coneSize(2),size(primaArray.activationDS,3));
    % Get center position of array
    centerOffset = [primaArray.height primaArray.width]/2;
    %     centerOffset = [2*primaArray.pixelWidth+primaArray.center(end,end,1)-primaArray.center(1,1,1) 2*primaArray.pixelWidth+primaArray.center(end,end,2)-primaArray.center(1,1,2)]/2;
    
    
    activationWindow = ceil([coneSize(1),coneSize(2)]/numberElectrodesX);

    % Set the spatial attenuation parameter
    shrinkFactor = 1;
    primaArray.spatialWeight = fspecial('Gaussian', round(1.25*activationWindow(1)), shrinkFactor*activationWindow(1)/2);
    primaArray.spatialWeight = primaArray.spatialWeight./max(primaArray.spatialWeight(:));
    
    % Compute full spatial scale of bp response due to each electrode
    szWeight = size(primaArray.spatialWeight);
    for ri = 1:size(primaArray.center,1)-0 % x-pos
        for ci = 1:size(primaArray.center,2)-0 % y-pos
            
            % Get center coords of each electrode
            centerCoords = scaleFactor*(centerOffset+[primaArray.center(ri,ci,1),primaArray.center(ri,ci,2)]);
            
            % Compute full spatial scale electrode coords
            primaArray.electrodeCoordsFull(ci,ri).x = centerCoords(1);
            primaArray.electrodeCoordsFull(ci,ri).y = centerCoords(2);
            primaArray.electrodeCoordsFull(ci,ri).rgb = primaArray.electrodeCoords(ci,ri).rgb;
            
            % Find indices of bipolar stimulated by this electrode
            rc = [ceil(-szWeight(1)/2+centerCoords(1)) : floor(szWeight(1)/2+centerCoords(1))];
            rcind1 = find(rc>0); rcind2 = find(rc<coneSize(1)); rcind =intersect(rcind1,rcind2); rc = rc(rcind);
            cc = [ceil(-szWeight(2)/2+centerCoords(2)) : floor(szWeight(2)/2+centerCoords(2))];
            ccind1 = find(cc>0); ccind2 = find(cc<coneSize(2)); ccind =intersect(ccind1,ccind2);  cc = cc(ccind);
            
            % Set spatial attenuation weight
            weightRS = primaArray.spatialWeight(rcind,ccind);
            allWeights{ci,ri} = weightRS';

            
            %             bpResponseCenterFull(rc,cc,:) = XW2RGBFormat(weightRS(:)*squeeze(bpResponseCenterTemp(ri,ci,:))',length(rc), length(cc));
            %             bpResponseSurroundFull(rc,cc,:) = XW2RGBFormat(weightRS(:)*squeeze(bpResponseSurroundTemp(ri,ci,:))',length(rc), length(cc));
            
            % Add each contribution of electrode to bipolar mosaic, weight
            % by spatial attenuation
            bpResponseCenterFull(rc,cc,:) = bpResponseCenterFull(rc,cc,:)+ XW2RGBFormat(weightRS(:)*squeeze(bpResponseCenterTemp(ri,ci,:))',length(rc), length(cc));
            bpResponseSurroundFull(rc,cc,:) = bpResponseSurroundFull(rc,cc,:) + XW2RGBFormat(weightRS(:)*squeeze(bpResponseSurroundTemp(ri,ci,:))',length(rc), length(cc));
            
        end
    end
    %     figure; imagesc(electrodeStimMask(:,:,10))
    
%     primaArray.visualizePhotocurrentAndBpMosaicResponses('bipolarActivation',...
%         allWeights, currentTemp, bpResponseCenterTemp, bpResponseCenterFull);
    
end

multFactor = [1 1 1 1];
for cellTypeInd = 1:4
    % Set into the bipolar mosaic
    bpL.mosaic{cellTypeInd}.set('responseCenter',multFactor(cellTypeInd)*bpResponseCenterFull);
    bpL.mosaic{cellTypeInd}.set('responseSurround',multFactor(cellTypeInd)*bpResponseSurroundFull);
end

% Set into the prima array
primaArray.bpMosaic = bpL;

end