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
eaSize = size(primaArray.center);
numberElectrodesX = eaSize(2);
numberElectrodesY = eaSize(1);
numberElectrodes = numberElectrodesX*numberElectrodesY;


%% Compute electrode activations from movieInput

% Get the full image/movie 
fullStimulus = permute(movieInput,[2 1 3]);

% Build the attenuation weighting from the center of the electrode;
activationWindow = ceil(size(movieInput,1)/numberElectrodesX);
electrodeAtten = fspecial('Gaussian', ceil(activationWindow), activationWindow/2);
electrodeAtten = electrodeAtten./max(electrodeAtten(:));
electrodeAttenTemporal = repmat(electrodeAtten,[1 1 size(fullStimulus,3)]);

sampleInterval = size(movieInput,1)/(numberElectrodesX);
% Elecrode coords in stimulus space
for xPos = 1:numberElectrodesX
    for yPos = 1:numberElectrodesY
        primaArray.electrodeCoords(xPos,yPos).x = (xPos-0.5)*sampleInterval;
        primaArray.electrodeCoords(xPos,yPos).y = (yPos-0.5)*sampleInterval;
        primaArray.electrodeCoords(xPos, yPos).rgb = [(xPos-1)/(numberElectrodesX-1)  0.2 (yPos-1)/(numberElectrodesY-1) ];
    end
end

% Find electrode activations by taking mean within spatial patch for each
% electrode over all time steps.
for xPos = 1:numberElectrodesX
    for yPos = 1:numberElectrodesY

        % Xcoords of window for stimulus
        imageCoordX1 = round(primaArray.electrodeCoords(xPos,yPos).x-activationWindow/2);
        imageCoordX2 = round(primaArray.electrodeCoords(xPos,yPos).x+activationWindow/2);
        
        % Ycoords of window for stimulus
        imageCoordY1 = round(primaArray.electrodeCoords(xPos,yPos).y-activationWindow/2);
        imageCoordY2 = round(primaArray.electrodeCoords(xPos,yPos).y+activationWindow/2);
        
        % Check to make sure we are not off of the edge of the stimulus
        imageCoordX1 = min([ floor(size(fullStimulus,2)-activationWindow/2)  max([1 imageCoordX1]) ]);
        imageCoordY1 = min([ floor(size(fullStimulus,1)-activationWindow/2)  max([1 imageCoordY1]) ]);
        imageCoordX2 = min([size(fullStimulus,2) imageCoordX2]);
        imageCoordY2 = min([size(fullStimulus,2) imageCoordY2]);
           
        % Pull out piece of stimulus and take mean
        electrodeStimulus = squeeze(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2,:,:));
        linearActivation(yPos, xPos,:) = mean(mean(electrodeStimulus,1),2);
        
        % primaArray.activation(xPos,yPos,:) = mean(RGB2XWFormat(electrodeStimulus));                
        
        % Implement the local electrode min([e1,e2]) nonlinearity
        sizeES = size(electrodeStimulus);
        
         if imageCoordX1 < (size(fullStimulus,2)-activationWindow/2)
            % electrodeAttenL = electrodeAttenTemporal(1:length(imageCoordY1:imageCoordY2), 1:length(imageCoordX1:floor(imageCoordX1+activationWindow/2)),:);
            % electrodeAttenR = electrodeAttenTemporal(1:length(imageCoordY1:imageCoordY2),end-length(floor(imageCoordX1+activationWindow/2)+1:imageCoordX2)+1:end,:);
            
            electrodeStimulusL = (fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:floor(imageCoordX1+activationWindow/2),:,:));
            electrodeStimulusR = (fullStimulus(imageCoordY1:imageCoordY2,floor(imageCoordX1+activationWindow/2)+1:imageCoordX2,:,:));
            % primaArray.activation(xPos,yPos,frame) = min([ mean(electrodeStimulus(:,1:floor(sizeES(2)/2))) mean(electrodeStimulus(:,ceil(sizeES(2)/2):sizeES(2)))]);
            % primaArray.activation(xPos,yPos,:) = min([mean(RGB2XWFormat(electrodeAttenL.*electrodeStimulusL)); mean(RGB2XWFormat(electrodeAttenR.*electrodeStimulusR))]);
            primaArray.activation(xPos,yPos,:) = min([mean(RGB2XWFormat(1.*electrodeStimulusL),1); mean(RGB2XWFormat(1.*electrodeStimulusR),1)]);
            
            %         imElec(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2) = (xPos*yPos)*ones(size(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2)));
        else
            
            % electrodeAttenL = electrodeAttenTemporal(1:length(imageCoordY1:imageCoordY2), 1:length(imageCoordX1:floor(imageCoordX2)),:);
            % electrodeAttenR = electrodeAttenTemporal(1:length(imageCoordY1:imageCoordY2),end-length(floor(imageCoordX1)+1:imageCoordX2)+1:end,:);
  
            electrodeStimulusL = (fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:floor(imageCoordX2),:,:));
            electrodeStimulusR = (fullStimulus(imageCoordY1:imageCoordY2,floor(imageCoordX1)+1:imageCoordX2,:,:)); 
%             szESR = size(((electrodeStimulusRsq)));
%             electrodeStimulusR = zeros(szESR(1),1,szESR(2));
%             electrodeStimulusR(1:szESR(1),:,1:szESR(2)) = electrodeStimulusRsq;
            % primaArray.activation(xPos,yPos,frame) = min([ mean(electrodeStimulus(:,1:floor(sizeES(2)/2))) mean(electrodeStimulus(:,ceil(sizeES(2)/2):sizeES(2)))]);
            % primaArray.activation(xPos,yPos,:) = min([mean(RGB2XWFormat(electrodeAttenL.*electrodeStimulusL)); mean(RGB2XWFormat(electrodeAttenR.*electrodeStimulusR))]);
            primaArray.activation(xPos,yPos,:) = min([mean(RGB2XWFormat(1.*electrodeStimulusL),1); mean(RGB2XWFormat(1.*electrodeStimulusR),1)]);
            
            %         imElec(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2) = (xPos*yPos)*ones(size(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2)));
        end
    end
end
% figure; imagesc(imElec);


%% Add X Hz spiking of stimulus
% primaArray.pulseFreq = pulseFreq;
% Right now the electrode sampling is at 0.008 s = 125 Hz
% Downsample if necessary
primaArray.pulseFreq
szAct = size(primaArray.activation);
primaArray.activationDS = zeros(szAct);
for iSample = 1:size(fullStimulus,3)
    if mod(iSample,100/primaArray.pulseFreq)==0
        primaArray.activationDS(:,:,iSample) = primaArray.activation(:,:,iSample);
        primaArray.activationDSoff(:,:,iSample) = primaArray.activation(:,:,iSample);
    end
end
ph=0;
% eaRS = reshape(primaArray.activation,[szAct(1)*szAct(2),szAct(3)]);
% eaDSRS = reshape(primaArray.activationDS,[szAct(1)*szAct(2),szAct(3)]);

%% Visualize stimulus and electrode activation
% primaArray.visualizeStimulusAndElectrodeActivation('electrodeActivation', fullStimulus, linearActivation, primaArray.activation, primaArray.activationDS, primaArray.activationDSoff)

end

