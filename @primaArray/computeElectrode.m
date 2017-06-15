function primaArray = computeElectrode(primaArray, movieInput)
%COMPUTEELECTRODE - a function of primaArray for computing the electrode activations.
% 
%   primaArray.computeElectrode(movieInput);
% 
% This function is called from primaArray.compute() and computes the
% activations of each electrode over time based on the movie input.
% 
% First, the size of the electrode array is determined, and then the
% relevant patch from movieInput is found for each electrode. There is a
% spatial falloff from each electrode so that the movie pixels closest to
% the center of the electrode have the largest effect on its activation
% level.
% 
% The local min(left_electrode,right_electrode) nonlinearity is also
% applied to the movieInput.
% 
% With the electrode activations over time determined, the signal is
% subsampled (if needed) to match the pulse frequency of the movie input.
% A negatively polarized version of the electrode activations is also
% captured in order to stimulate off bipolar cells.
% 
% 
% 5/2017 JRG (c) isetbio team
    
%% Compute the size of the electrode array 

numberElectrodesX = floor(primaArray.width/primaArray.pixelWidth)+1;
numberElectrodesY = floor(primaArray.width/primaArray.pixelWidth)+1;
numberElectrodes = numberElectrodesX*numberElectrodesY;

activationWindow = ceil(size(movieInput,1)/numberElectrodesX);

%% Compute electrode activations from movieInput

% Get the full image/movie from the identity outersegment
fullStimulus = movieInput;

% Build the attenuation weighting from the center of the electrode;
electrodeAtten = fspecial('Gaussian', round(activationWindow(1)), activationWindow(1)/2);
electrodeAtten = electrodeAtten./max(electrodeAtten(:));
electrodeAttenTemporal = repmat(electrodeAtten,[1 1 size(fullStimulus,3)]);

% Find electrode activations by taking mean within spatial patch for each
% electrode over all time steps.
for xPos = 1:numberElectrodesX
    for yPos = 1:numberElectrodesY

        % Xcoords of window for stimulus
        imageCoordX1 = (activationWindow)*(xPos-1)+1;
        imageCoordX2 = (activationWindow)*(xPos);
        
        % Ycoords of window for stimulus
        imageCoordY1 = (activationWindow)*(yPos-1)+1;
        imageCoordY2 = (activationWindow)*(yPos);
        
        % Check to make sure we are not off of the edge of the stimulus
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
            % electrodeAttenL = electrodeAttenTemporal(1:length(imageCoordY1:imageCoordY2), 1:length(imageCoordX1:floor(imageCoordX1+activationWindow/2)),:);
            % electrodeAttenR = electrodeAttenTemporal(1:length(imageCoordY1:imageCoordY2),end-length(floor(imageCoordX1+activationWindow/2)+1:imageCoordX2)+1:end,:);
            
            electrodeStimulusL = squeeze(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:floor(imageCoordX1+activationWindow/2),:,:));
            electrodeStimulusR = squeeze(fullStimulus(imageCoordY1:imageCoordY2,floor(imageCoordX1+activationWindow/2)+1:imageCoordX2,:,:));
            % primaArray.activation(xPos,yPos,frame) = min([ mean(electrodeStimulus(:,1:floor(sizeES(2)/2))) mean(electrodeStimulus(:,ceil(sizeES(2)/2):sizeES(2)))]);
            % primaArray.activation(xPos,yPos,:) = min([mean(RGB2XWFormat(electrodeAttenL.*electrodeStimulusL)); mean(RGB2XWFormat(electrodeAttenR.*electrodeStimulusR))]);
            primaArray.activation(xPos,yPos,:) = min([mean(RGB2XWFormat(1.*electrodeStimulusL)); mean(RGB2XWFormat(1.*electrodeStimulusR))]);
            
            %         imElec(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2) = (xPos*yPos)*ones(size(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2)));
        else
            
            % electrodeAttenL = electrodeAttenTemporal(1:length(imageCoordY1:imageCoordY2), 1:length(imageCoordX1:floor(imageCoordX2)),:);
            % electrodeAttenR = electrodeAttenTemporal(1:length(imageCoordY1:imageCoordY2),end-length(floor(imageCoordX1)+1:imageCoordX2)+1:end,:);
  
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
% Right now the electrode sampling is at 0.008 s = 125 Hz
% Downsample if necessary

szAct = size(primaArray.activation);
primaArray.activationDS = zeros(szAct);
for iSample = 1:size(fullStimulus,3)
    if mod(iSample,100/primaArray.pulseFreq)==0
        primaArray.activationDS(:,:,iSample) = primaArray.activation(:,:,iSample);
        primaArray.activationDSoff(:,:,iSample) = 1-primaArray.activation(:,:,iSample);
    end
end

% eaRS = reshape(primaArray.activation,[szAct(1)*szAct(2),szAct(3)]);
% eaDSRS = reshape(primaArray.activationDS,[szAct(1)*szAct(2),szAct(3)]);