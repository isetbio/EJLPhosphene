% 
% A first-pass attempt at simulating the activations of retinal
% ganglion cells with an array of stimulating electrodes. The stimulus
% image is defined, the electrode array is generated, the electrode
% activations are computed, the RGC mosaics are generated, the RGC mosaic
% responses are computed and the stimulus is inferred from the RGC mosaic
% resposnes using simple linear summation of the STAs.
% 
% Outline of computation:
% 1. Load image/movie
% - Moving bar stimulus sweeping left to right
% 
% 2. Outer segment representation
% - The osIdentity object does no computation and stores the stimulus RGB
% data.
%  Check iestimulusbar for size change
% 
% 3. Build electrode array
% - The electrode array is simulated as a square array of a certain size
% with Gaussian weights on their activations. The activation of each
% electrode is the mean of the CHANGE TO HEX
% 
% 4. Compute electrode activations from image/movie
% 5. Build RGC array
% 6. Calculate RGC input - only one, change to input from multiple
% 7. Build RGC activation functions
% 8. Compute RGC activations/spikes
% 9. Invert representation to form image/movie
% 
% 3/2016 JRG (c) isetbio
% 
% 1.5 orders of magnitude variance in threshold of RGC activation
% 
% Add 5 Hz stimulation pulses




%% Initialize
% clear;
% ieInit;

tic

bwL = 20%[ 32 50]
freqL = 5%[5 8 2 12]

%% Parameters to alter
clear electrodeArray
% Electrode size
% Set the size of implant pixels
electrodeArray.width = 1*15e-6; % meters
% electrodeArray.width = 140e-6; % meters
 
% Retinal patch eccentricity
patchEccentricity = 4; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.6*2/3;

% Stimulus length
nSteps = 100;

% Activation curve

% Spatial activation of electrode

% Electrode PWM 

percentDead = 0;
%% Load image
clear params
% One frame of a moving bar stimulus
% Set parameters for size
% params.nSteps = nSteps;
params.barWidth = 20;
params.row = 96;
params.col = 96;
params.fov = fov;
% params.freq = freqL; % Hz grating frequency
% % params.vfov = 0.7;
params.startFrames = 0;
params.endFrames = 0;
movingBar = ieStimulusBarSensor(params);
iStim = movingBar;
tuningWoffElec = 1;
tuningWoffHealthy = 1;

pulseFreq = 25; % Hz, electrode pulse frequency

contrastHealthy = 1;
contrastElectrode = 1;
%%% Grating subunit stimulus

% params.barWidth = bwL;
% iStim = ieStimulusGratingSubunit(params);

% iStim = iStimC;
absorptions = iStim.absorptions;
% movingBar = iStim;

nSteps = size(iStim.sceneRGB,3);
params.nSteps = nSteps;

size(iStim.sceneRGB)
%% Show raw stimulus for osIdentity
% figure;
% for frame1 = 1:size(movingBar.sceneRGB,3)
%     imagesc(squeeze(movingBar.sceneRGB(:,:,frame1,:)));
%     colormap gray; drawnow;
% end
% close;

%% Outer segment calculation
% There is no simulated outer segment, this identity outer segment acts as
% a pass-through holder of the stimulus intensity information.

% Input = RGB
os = osCreate('displayrgb');

sceneSize = sceneGet(movingBar.scene,'size');
retinalPatchWidth = sensorGet(movingBar.absorptions,'width','m');
% retinalPatchWidth = sceneGet(movingBar.scene,'width');

% retinalPatchHeight = sensorGet(movingBar.absorptions,'height','m');
retinalPatchHeight = (sceneSize(1)/sceneSize(2))*retinalPatchWidth;

% % % coneSpacing = scene.wAngular*300
% coneSpacing = sensorGet(sensor,'dimension','um');
os = osSet(os, 'patchSize', retinalPatchWidth);

timeStep = sensorGet(movingBar.absorptions,'time interval','sec');
os = osSet(os, 'timeStep', timeStep);

% movingBar.sceneRGB = (contrastElectrode)*(movingBar.sceneRGB - 0.5)+0.5;
sceneRGB_pros = (contrastElectrode)*(movingBar.sceneRGB - 0.5)+0.5;
os = osSet(os, 'rgbData', sceneRGB_pros);

sceneRGB_Healthy = (contrastHealthy)*(movingBar.sceneRGB - 0.5)+0.5;
osHealthy = os;
osHealthy = osSet(osHealthy, 'rgbData', sceneRGB_Healthy);

% os = osCompute(absorptions);

% % Plot the photocurrent for a pixel
% osPlot(os,absorptions);

retinalPatchSize = osGet(os,'size');
numberElectrodesX = floor(retinalPatchWidth/electrodeArray.width)+4;
numberElectrodesY = floor(retinalPatchWidth/electrodeArray.width)+0;
numberElectrodes = numberElectrodesX*numberElectrodesY;
%% Build electrode array
% Define the electrode array structure/object

% Size stores the size of the array of electrodes
electrodeArray.size = [numberElectrodesX numberElectrodesY];

% Build the matrix of center coordinates for each electrode
% electrodeArray.center(xPos,yPos,:) = [xCoord yCoord];
x0 = -(numberElectrodesX-1)*electrodeArray.width/2;
y0 = -(numberElectrodesY-1)*electrodeArray.width/2;
for xPos = 1:numberElectrodesX
    for yPos = 1:numberElectrodesY
        % electrodeArray.center(xPos,yPos,:) = [x0+(electrodeArray.width/2)*(xPos-1) + electrodeArray.width, y0+(electrodeArray.width/2)*(yPos-1) + electrodeArray.width];
        electrodeArray.center(xPos,yPos,:) = [x0+(electrodeArray.width/1)*(xPos-1) + 0, y0+(electrodeArray.width/1)*(yPos-1) + 0 + (mod(xPos,2)-.5)*(electrodeArray.width/2)];
    end
end
% xe = electrodeArray.center(:,:,1); ye = electrodeArray.center(:,:,2);
% figure; scatter(xe(:),ye(:));

th = (0:1/6:1)'*2*pi;
xh = electrodeArray.width/2*cos(th);
yh = electrodeArray.width/2*sin(th);
% figure
% fill(x,y,'r')
% axis square

% % Plot electrode array
eaSize = size(electrodeArray.center);
figure;
hold on;
for i = 1:eaSize(1)
    for j = 1:eaSize(2)
%         scatter(electrodeArray.center(i,j,1),electrodeArray.center(i,j,2));
        plot(xh+electrodeArray.center(i,j,1),yh+electrodeArray.center(i,j,2),'r')
    end
end
axis equal
xlabel('Distance (m)'); ylabel('Distance (m)');
set(gca,'fontsize',14);

% Build the current stimulation activation window
% Gaussian activation from center of electrode
activationWindow = floor(retinalPatchSize(2)/numberElectrodesX);
electrodeArray.spatialWeight = fspecial('Gaussian', activationWindow, activationWindow/8);

% Visualize Gaussian activation
% figure; imagesc(electrodeArray.spatialWeight); 
% figure; surf(electrodeArray.spatialWeight); 
% xlabel(sprintf('Distance (\\mum)')); ylabel(sprintf('Distance (\\mum)'));
% title('Gaussian Activation for a Single Electrode'); set(gca,'fontsize',16);
%% Compute electrode activations from image

% Get the full image/movie from the identity outersegment
fullStimulus = osGet(os,'rgbData');

% Find electrode activations by taking mean within window
for frame = 1:params.nSteps
    for xPos = 1:numberElectrodesX
        for yPos = 1:numberElectrodesY
            % Xcoords of window for stimulus
            imageCoordX1 = (activationWindow)*(xPos-1)+1;
            imageCoordX2 = (activationWindow)*(xPos);
            
            % Ycoords of window for stimulus
            imageCoordY1 = (activationWindow)*(yPos-1)+1;
            imageCoordY2 = (activationWindow)*(yPos);
            
            if imageCoordX2 > size(fullStimulus,2); imageCoordY2 = size(fullStimulus,2); end;
            if imageCoordY2 > size(fullStimulus,1); imageCoordY2 = size(fullStimulus,1); end;
            % Pull out piece of stimulus and take mean
            electrodeStimulus = squeeze(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2,frame,:));
            electrodeArray.activation(xPos,yPos,frame) = mean(electrodeStimulus(:));
            
%             sizeES = size(electrodeStimulus);
%             electrodeArray.activation(xPos,yPos,frame) = min([ mean(electrodeStimulus(:,1:floor(sizeES(2)/2))) mean(electrodeStimulus(:,ceil(sizeES(2)/2):sizeES(2)))]);

            % imagesc(electrodeStimulus); title(sprintf('%2.2f',mean(electrodeStimulus(:))));
        end
    end
end

% Show plot

% Apply Gaussian

%% Just electrode activation
% With < 20 uC/electrode there was no perception, Zrenner paper

% % Col activation
% electrodeArray.activation(:) = 0;
% electrodeArray.activation(4,:,4:8) = 0.1;
% electrodeArray.activation(4,:,24:28) = 0.5;
% electrodeArray.activation(4,:,44:48) = 0.75;
% electrodeArray.activation(4,:,64:68) = 1;

% name_str = 'col_electrode_10fps_lowOFF.mp4';
% electrodeArray.activation(:) = 0;
% frame = 20;
%     for xPos = 2:2:numberElectrodesX-2
% %         for yPos = 2:2:numberElectrodesY
%             frame = frame + 20;
%             electrodeArray.activation(xPos,:,frame:frame+4) = 1;
%             
% %         end
%     end
%     electrodeArray.activation(xPos,yPos,frame+20) = 0;
% params.nSteps = frame+20;

%%% Single electrode activation
% name_str = 'single_electrode_10fps_lowOFF.mp4';
% electrodeArray.activation(:) = 0;
% frame = 20;
%     for xPos = 2:2:numberElectrodesX
%         for yPos = 2:2:numberElectrodesY
%             frame = frame + 20;
%             electrodeArray.activation(xPos,yPos,frame:frame+4) = 1;
%             
%         end
%     end
%     electrodeArray.activation(xPos,yPos,frame+20) = 0;
% params.nSteps = frame+20;


% % % U activation
% name_str = 'U_10fps.mp4';
% electrodeArray.activation(:) = 0;
% frame = 20;
% for rep = 1:3
%     electrodeArray.activation(3,1:4,frame+1:frame+4) = 1;
%     electrodeArray.activation(6,1:4,frame+1:frame+4) = 1;
%     electrodeArray.activation(3:6,4,frame+1:frame+4) = 1;
%     frame = frame+30;
% end


%%%%%%
% electrodeArray.activation(4,:,4:8) = 0.1;
% electrodeArray.activation(4,:,24:28) = 0.5;
% electrodeArray.activation(4,:,44:48) = 0.75;
% electrodeArray.activation(4,:,64:68) = 1;
%% Add 5 Hz spiking of stimulus

% Right now the electrode sampling is at 0.01 s = 100 Hz
% Downsample to get 5 Hz
szAct = size(electrodeArray.activation);
electrodeArray.activationDS = zeros(szAct);
for iSample = 1:szAct(3)
    if mod(iSample,100/pulseFreq)==0
    electrodeArray.activationDS(:,:,iSample) = electrodeArray.activation(:,:,iSample);
    electrodeArray.activationDSoff(:,:,iSample) = 1-electrodeArray.activation(:,:,iSample);
    end
end

eaRS = reshape(electrodeArray.activation,[szAct(1)*szAct(2),szAct(3)]);
eaDSRS = reshape(electrodeArray.activationDS,[szAct(1)*szAct(2),szAct(3)]);
figure; 
% plot(eaRS'); 
hold on; 
plot(eaDSRS');

%% Build RGC array

clear paramsIR innerRetina
paramsIR.name    = 'Macaque inner retina pixium 1'; % This instance
paramsIR.eyeSide   = 'left';   % Which eye
paramsIR.eyeRadius = 3;        % Radius in mm
paramsIR.eyeAngle  = 90;       % Polar angle in degrees

% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/WNstim_response_OffParasol_64_grating_june10.mat');
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/WNstim_response_OffParasol_RGC.mat')
% innerRetina2 = innerRetina; clear innerRetina;
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/WNstim_response_OnParasol_RGC.mat')



filenameRGC = ['/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/WNstim_response_OffMidget_RGC.mat'];
load(filenameRGC);
innerRetina3 = innerRetina;

filenameRGC = ['/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/WNstim_response_OnMidget_RGC.mat'];
load(filenameRGC);
innerRetina4 = innerRetina;

rdt = RdtClient('isetbio');
rdt.crp('/resources/data/rgc/pixium/');
data = rdt.readArtifact('WNstim_response_OffParasol_RGC', 'type', 'mat');
innerRetina2 = data.innerRetina; % clear innerRetina;
data = rdt.readArtifact('WNstim_response_OnParasol_RGC', 'type', 'mat');
innerRetina = data.innerRetina; 

% filenameRGC = ['C:\Users\James\Documents\GitHub\May26_onBig2\WNstim_response_OnParasol_RGC.mat'];

% filenameRGC = ['C:\Users\James\Documents\GitHub\May26_offBig1\WNstim_response_OffParasol_RGC.mat'];


% model   = 'LNP';    % Computational model
% innerRetina = irCreate(os,paramsIR);
% innerRetina = rgcMosaicCreate(innerRetina,'type','onMidget','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','offMidget','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','onParasol','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','offParasol','model',model);
% % innerRetina = rgcMosaicCreate(innerRetina,'type','sbc','model',model);


% % % % % Plot
% irPlot(innerRetina,'mosaic');
% % % figure;
% hold on;
% for spInd = 1:length(innerRetina.mosaic)
% for i = 3:eaSize(1)-2
%     for j = 1:eaSize(2)
%         subplot(floor(sqrt(length(innerRetina.mosaic))),ceil(sqrt(length(innerRetina.mosaic))),spInd); 
%         hold on;
% %         scatter(electrodeArray.center(i,j,1),electrodeArray.center(i,j,2));
%         plot(xh+electrodeArray.center(i,j,1),yh+electrodeArray.center(i,j,2),'r','linewidth',3)
%     end
% end
% end
% 
metersPerPixel = retinalPatchWidth/retinalPatchSize(2);
% 
% % innerRetina = irCompute(innerRetina, os);
% 
% % Plot the electrode array over spatial receptive fields
% % eaSize = size(electrodeArray.center);
% % hold on;
% % for i = 1:eaSize(1)
% %     for j = 1:eaSize(2)
% %         scatter(1e7*electrodeArray.center(i,j,1),1e7*electrodeArray.center(i,j,2));
% %     end
% % end
% mosaicInd = 1;
% 
electrodeCenter = reshape(electrodeArray.center,[numberElectrodes,2]);

rp = metersPerPixel*vertcat(innerRetina.mosaic{1}.cellLocation{:});
% figure; scatter(rp(:,1),rp(:,2))
% hold on; scatter(electrodeCenter(:,1),electrodeCenter(:,2),'xr')

%% Tile mosaics for reconstruction

electrodePos0 = electrodeCenter(1,:);
electrodePos1 = electrodeCenter(end,:);

rgcPos0 = rp(1,:); rgcPos1 = rp(end,:);

electrodeDistRows = (electrodePos1(1) - electrodePos0(1));
electrodeDistCols = (electrodePos1(2) - electrodePos0(2));
rgcDistRows = (rgcPos1(1) - rgcPos0(1));
rgcDistCols = (rgcPos1(2) - rgcPos0(2));

nTileRows = round( electrodeDistRows / rgcDistRows);
nTileCols = round( electrodeDistCols / rgcDistCols);

mosaicCtr = 0;
for iTileRows = 1:nTileRows
    for iTileCols = 1:nTileCols
        mosaicCtr = mosaicCtr+1;
        innerRetina.mosaic{mosaicCtr} = innerRetina.mosaic{1};
        innerRetina2.mosaic{mosaicCtr} = innerRetina2.mosaic{1};
        
        innerRetina3.mosaic{mosaicCtr} = innerRetina3.mosaic{1};
        
        innerRetina4.mosaic{mosaicCtr} = innerRetina4.mosaic{1};
        mosaicOffset(iTileRows, iTileCols, 1) = electrodePos0(1) + (rgcDistRows)*(iTileRows-1) + (rgcDistRows/2);
        mosaicOffset(iTileRows, iTileCols, 2) = electrodePos0(2) + (rgcDistCols)*(iTileCols-1) +(rgcDistCols/2);
 
    end
end

size(mosaicOffset)
[nTileRows*nTileCols 2]
% mosaicOffsetPlot = reshape(mosaicOffset, [nTileRows*nTileCols 2]);
% hold on; scatter(mosaicOffsetPlot(:,1), mosaicOffsetPlot(:,2),'go','filled');

%% Calculate RGC input
% Weight electrode activation by Gaussian as a function of distance between
% centers
params.nSteps = params.nSteps-5;
% On midget
offFlag = 0;
innerRetinaInput = irActivationFromElectrode(innerRetina, electrodeArray, retinalPatchWidth, metersPerPixel, nTileRows, nTileCols, mosaicOffset, params, offFlag);

% Off parasol
offFlag = 1;
innerRetinaInput2 = irActivationFromElectrode(innerRetina2, electrodeArray, retinalPatchWidth, metersPerPixel, nTileRows, nTileCols, mosaicOffset, params, offFlag);


% Off midget
offFlag = 0;
innerRetinaInput3 = irActivationFromElectrode(innerRetina3, electrodeArray, retinalPatchWidth, metersPerPixel, nTileRows, nTileCols, mosaicOffset, params, offFlag);


% On midget
offFlag = 0;
innerRetinaInput4 = irActivationFromElectrode(innerRetina4, electrodeArray, retinalPatchWidth, metersPerPixel, nTileRows, nTileCols, mosaicOffset, params, offFlag);

%% Build RGC activation functions

innerRetinaThreshold = irGetThreshold(innerRetina);

innerRetinaThreshold2 = irGetThreshold(innerRetina2);

innerRetinaThreshold3 = irGetThreshold(innerRetina3);

innerRetinaThreshold4 = irGetThreshold(innerRetina4);

%% Compute RGC activations

innerRetina = irGetLinearRespElectrode(innerRetina, 10*innerRetinaInput, innerRetinaThreshold, params);
innerRetina2 = irGetLinearRespElectrode(innerRetina2, 10*innerRetinaInput2, innerRetinaThreshold2, params);
innerRetina3 = irGetLinearRespElectrode(innerRetina3, 10*innerRetinaInput3, innerRetinaThreshold3, params);
innerRetina4 = irGetLinearRespElectrode(innerRetina4, 10*innerRetinaInput4, innerRetinaThreshold4, params);


% irPlot(innerRetina,'mosaic');
% % Visualize thresholds
% figure; hold on;
% for ii = 1:xc
%     for ji = 1:yc
%         xp = -.5:.01:1;
%         plot((xp),1./(1+exp(-thr*(i0all{mosaicInd}(ii,ji)+(xp)))));
%     end
% end

% xpn = -.5:.01:4;
% figure;
% plot((xpn),1./(1+exp(-2*(-0.6 + (xpn)))));
% 
% % plot(eaDSRS');
% % figure; plot(1./(1+exp(-2*(-0.6 + (eaDSRS)))));
% 
% figure; plot(1./(1+exp(-2*(-0.6 + (innerRetina.mosaic{1}.responseLinear{1,1})))));

%% Compute RGC spiking
numberTrials = 1;
for tr = 1:numberTrials
    innerRetina = irComputeSpikes(innerRetina,'coupling',false);
    innerRetina2 = irComputeSpikes(innerRetina2,'coupling',false);
    innerRetina3 = irComputeSpikes(innerRetina3,'coupling',false);
    innerRetina4 = irComputeSpikes(innerRetina4,'coupling',false);
end
clear innerRetina1
innerRetina1 = innerRetina;
% irPlot(innerRetina, 'linear');
%% Invert representation to form image/movie
clear stimulusReconstruction

m1 = min(min(mosaicOffset(:,:,1)./metersPerPixel));
m2 = min(min(mosaicOffset(:,:,2)./metersPerPixel));

offsetStitch = mosaicOffset./metersPerPixel;
offsetStitch(:,:,1) = offsetStitch(:,:,1) - m1;
offsetStitch(:,:,2) = offsetStitch(:,:,2) - m2;
% [stimulusReconstruction, stimulusReconstructionOnes, paramsRec] = irReconstructStitch(innerRetina, 'tuningWoff', 1,'mosaicOffset',offsetStitch);

% figure; ieMovie(stimulusReconstruction(1:100,1:100,:));

%%
clear movrecons_on_off movieStitch
[movrecons_on_off, movrecons_on_off_dropout] = irOptimalRecon(innerRetina1, innerRetina2, innerRetina3, innerRetina4, percentDead);


%% Stitch together recon movie from different mosaics

% electrodePos0 = electrodeCenter(1,:);
% electrodePos1 = electrodeCenter(end,:);
% 
% rgcPos0 = rp(1,:); rgcPos1 = rp(end,:);
% 
% electrodeDistRows = (electrodePos1(1) - electrodePos0(1));
% electrodeDistCols = (electrodePos1(2) - electrodePos0(2));
% rgcDistRows = (rgcPos1(1) - rgcPos0(1));
% rgcDistCols = (rgcPos1(2) - rgcPos0(2));
% 
% nTileRows = round( electrodeDistRows / rgcDistRows);
% nTileCols = round( electrodeDistCols / rgcDistCols);

movieSize = size(movrecons_on_off{1});
% movieStitch = zeros(nTileRows*96,nTileCols*96,movieSize(3));
% movieStitchDropout = zeros(nTileRows*96,nTileCols*96,movieSize(3));

% nCut = 19;
nCut = 1;
movieStitch = zeros(nTileRows*96,nTileCols*96-nTileCols*nCut,movieSize(3));
movieStitchDropout = zeros(nTileRows*96,nTileCols*96-nTileCols*nCut,movieSize(3));

mosaicCtr = 0;
for iTileRows = 1:nTileRows
    for iTileCols = 1:nTileCols
        mosaicCtr = mosaicCtr+1;
%         innerRetina.mosaic{mosaicCtr} = innerRetina.mosaic{1};
%         innerRetina2.mosaic{mosaicCtr} = innerRetina2.mosaic{1};
%         mosaicOffset(iTileRows, iTileCols, 1) = electrodePos0(1) + (rgcDistRows)*(iTileRows-1) + (rgcDistRows/2);
%         mosaicOffset(iTileRows, iTileCols, 2) = electrodePos0(2) + (rgcDistCols)*(iTileCols-1) +(rgcDistCols/2);
        
%         movieStitch(1+(iTileRows-1)*96:96+(iTileRows-1)*96,1+(iTileCols-1)*96:96+(iTileCols-1)*96,:) = movrecons_on_off{mosaicCtr};
%         movieStitchDropout(1+(iTileRows-1)*96:96+(iTileRows-1)*96,1+(iTileCols-1)*96:96+(iTileCols-1)*96,:) = movrecons_on_off_dropout{mosaicCtr};
        
        
        movieStitch(1+(iTileRows-1)*96:96+(iTileRows-1)*96,1+(iTileCols-1)*96  - (iTileCols-1)*(nCut-1):96+(iTileCols-1)*96 - iTileCols*(nCut-1),:) = movrecons_on_off{mosaicCtr}(:,1:end-(nCut-1),:);
        movieStitchDropout(1+(iTileRows-1)*96:96+(iTileRows-1)*96 ,1+(iTileCols-1)*96- (iTileCols-1)*(nCut-1):96+(iTileCols-1)*96 - iTileCols*(nCut-1),:) = movrecons_on_off_dropout{mosaicCtr}(:,1:end-(nCut-1),:);
    end
end

figure; ieMovie(movieStitch);
% figure; ieMovie(movieStitchHealthy);

%% Build RGC array for healthy retina

clear paramsIR innerRetinaHealthy
% paramsIR.name    = 'Macaque inner retina 1'; % This instance
% paramsIR.eyeSide   = 'left';   % Which eye
% paramsIR.eyeRadius = 3;        % Radius in mm
% paramsIR.eyeAngle  = 90;       % Polar angle in degrees
% 
% model   = 'LNP';    % Computational model
% innerRetinaHealthy = irCreate(osHealthy,paramsIR);
% % innerRetinaHealthy = rgcMosaicCreate(innerRetinaHealthy,'type','onMidget','model',model);
% % innerRetinaHealthy = rgcMosaicCreate(innerRetinaHealthy,'type','offMidget','model',model);
% % innerRetinaHealthy = rgcMosaicCreate(innerRetinaHealthy,'type','onParasol','model',model);
% innerRetinaHealthy = rgcMosaicCreate(innerRetinaHealthy,'type','offParasol','model',model);


filenameRGC = ['/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/WNstim_response_OffMidget_RGC.mat'];
load(filenameRGC);
innerRetinaHealthy3 = innerRetina;

filenameRGC = ['/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/WNstim_response_OnMidget_RGC.mat'];
load(filenameRGC);
innerRetinaHealthy4 = innerRetina;

rdt = RdtClient('isetbio');
rdt.crp('/resources/data/rgc/pixium/');
data = rdt.readArtifact('WNstim_response_OffParasol_RGC', 'type', 'mat');
innerRetinaHealthy2 = data.innerRetina; % clear innerRetina;
data = rdt.readArtifact('WNstim_response_OnParasol_RGC', 'type', 'mat');
innerRetinaHealthy = data.innerRetina; 

%%
for mosaicInd = 2:length(innerRetina2.mosaic)
    innerRetinaHealthy.mosaic{mosaicInd} = innerRetinaHealthy.mosaic{1};
    szMos = innerRetinaHealthy.mosaic{mosaicInd}.get('mosaic size');
    
    [mosaicSubRow, mosaicSubCol] = ind2sub([nTileRows nTileCols],mosaicInd);
    newLocs = cell(szMos(1),szMos(2));
    for xc = 1:szMos(1)
        for yc = 1:szMos(2)
            newLocs{xc,yc} = ...
                innerRetinaHealthy.mosaic{mosaicInd}.cellLocation{xc,yc} + offsetStitch(mosaicSubRow,mosaicSubCol);
        end
    end
    innerRetinaHealthy.mosaic{mosaicInd}.set('cellLocation',newLocs);

end

for mosaicInd = 2:length(innerRetina2.mosaic)
    innerRetinaHealthy2.mosaic{mosaicInd} = innerRetinaHealthy2.mosaic{1};
    szMos = innerRetinaHealthy2.mosaic{mosaicInd}.get('mosaic size');
    
    [mosaicSubRow, mosaicSubCol] = ind2sub([nTileRows nTileCols],mosaicInd);
    newLocs = cell(szMos(1),szMos(2));
    for xc = 1:szMos(1)
        for yc = 1:szMos(2)
            newLocs{xc,yc} = ...
                innerRetinaHealthy2.mosaic{mosaicInd}.cellLocation{xc,yc} + offsetStitch(mosaicSubRow,mosaicSubCol);
        end
    end
    innerRetinaHealthy2.mosaic{mosaicInd}.set('cellLocation',newLocs);

end

for mosaicInd = 2:length(innerRetina2.mosaic)
    innerRetinaHealthy3.mosaic{mosaicInd} = innerRetinaHealthy3.mosaic{1};
    szMos = innerRetinaHealthy3.mosaic{mosaicInd}.get('mosaic size');
    
    [mosaicSubRow, mosaicSubCol] = ind2sub([nTileRows nTileCols],mosaicInd);
    newLocs = cell(szMos(1),szMos(2));
    for xc = 1:szMos(1)
        for yc = 1:szMos(2)
            newLocs{xc,yc} = ...
                innerRetinaHealthy3.mosaic{mosaicInd}.cellLocation{xc,yc} + offsetStitch(mosaicSubRow,mosaicSubCol);
        end
    end
    innerRetinaHealthy3.mosaic{mosaicInd}.set('cellLocation',newLocs);

end

for mosaicInd = 2:length(innerRetina2.mosaic)
    innerRetinaHealthy4.mosaic{mosaicInd} = innerRetinaHealthy4.mosaic{1};
    szMos = innerRetinaHealthy4.mosaic{mosaicInd}.get('mosaic size');
    
    [mosaicSubRow, mosaicSubCol] = ind2sub([nTileRows nTileCols],mosaicInd);
    newLocs = cell(szMos(1),szMos(2));
    for xc = 1:szMos(1)
        for yc = 1:szMos(2)
            newLocs{xc,yc} = ...
                innerRetinaHealthy4.mosaic{mosaicInd}.cellLocation{xc,yc} + offsetStitch(mosaicSubRow,mosaicSubCol);
        end
    end
    innerRetinaHealthy4.mosaic{mosaicInd}.set('cellLocation',newLocs);

end

%%
innerRetinaHealthy = irCompute(innerRetinaHealthy,osHealthy,'coupling',false);
innerRetinaHealthy2 = irCompute(innerRetinaHealthy2,osHealthy,'coupling',false);
innerRetinaHealthy3 = irCompute(innerRetinaHealthy3,osHealthy,'coupling',false);
innerRetinaHealthy4 = irCompute(innerRetinaHealthy4,osHealthy,'coupling',false);

%%

% clear movrecons_on_off movieStitch
[movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalRecon(innerRetinaHealthy, innerRetinaHealthy2, innerRetinaHealthy3, innerRetinaHealthy4, percentDead);

%% Movie stitch healthy
movieSize = size(movrecons_on_offHealthy{1});
% movieStitchHealthy = zeros(nTileRows*96,nTileCols*96,movieSize(3));
% movieStitchHealthyDropout = zeros(nTileRows*96,nTileCols*96,movieSize(3));


% nCut = 1;
movieStitchHealthy = zeros(nTileRows*96,nTileCols*96-nTileCols*nCut,movieSize(3));
movieStitchHealthyDropout = zeros(nTileRows*96,nTileCols*96-nTileCols*nCut,movieSize(3));

mosaicCtr = 0;
for iTileRows = 1:nTileRows
    for iTileCols = 1:nTileCols
        mosaicCtr = mosaicCtr+1;
%         innerRetina.mosaic{mosaicCtr} = innerRetina.mosaic{1};
%         innerRetina2.mosaic{mosaicCtr} = innerRetina2.mosaic{1};
%         mosaicOffset(iTileRows, iTileCols, 1) = electrodePos0(1) + (rgcDistRows)*(iTileRows-1) + (rgcDistRows/2);
%         mosaicOffset(iTileRows, iTileCols, 2) = electrodePos0(2) + (rgcDistCols)*(iTileCols-1) +(rgcDistCols/2);
        
%         movieStitchHealthy(1+(iTileRows-1)*96:96+(iTileRows-1)*96,1+(iTileCols-1)*96:96+(iTileCols-1)*96,:) = movrecons_on_offHealthy{mosaicCtr};
%         movieStitchHealthyDropout(1+(iTileRows-1)*96:96+(iTileRows-1)*96,1+(iTileCols-1)*96:96+(iTileCols-1)*96,:) = movrecons_on_offHealthy_dropout{mosaicCtr};
        
        movieStitchHealthy(1+(iTileRows-1)*96:96+(iTileRows-1)*96,1+(iTileCols-1)*96  - (iTileCols-1)*(nCut-1):96+(iTileCols-1)*96 - iTileCols*(nCut-1),:) = movrecons_on_offHealthy{mosaicCtr}(:,1:end-(nCut-1),:);
        movieStitchHealthyDropout(1+(iTileRows-1)*96:96+(iTileRows-1)*96 ,1+(iTileCols-1)*96- (iTileCols-1)*(nCut-1):96+(iTileCols-1)*96 - iTileCols*(nCut-1),:) = movrecons_on_offHealthy_dropout{mosaicCtr}(:,1:end-(nCut-1),:);
 
    end
end

figure; ieMovie(movieStitchHealthy);

toc
%%
% name_str = ['gratingH_20Hz_width_' num2str(params.barWidth) '_freq_' num2str(freqL) '_onM_8_hz_ON_IMS1_' num2str(cputime*100) '.mp4'];
name_str = ['prosthesis_combined_bar_' num2str(cputime*100) '.mp4'];
path_str = '/Users/james/Documents/MATLAB/isetbio misc/pixium_videos/meeting_may27/';
vObj = VideoWriter([path_str name_str],'MPEG-4');
vObj.FrameRate = 10;
vObj.Quality = 100;
open(vObj);

sizeScene = size(movingBar.sceneRGB(:,:,1,:));


caxispros = [1*min(movieStitch(:)) 1*max(movieStitch(:))];
caxishealthy = [1*min(movieStitchHealthy(:)) 1*max(movieStitchHealthy(:))];
% caxispros = [1*min(movieStitchDropout(:)) 1*max(movieStitchDropout(:))];
% caxishealthy = [1*min(movieStitchHealthyDropout(:)) 1*max(movieStitchHealthyDropout(:))];

% Play the movie with the stimulus
for loopv = 1%:10
h1=figure; 
% set(gcf,'position',[160 60 1070 740]);
set(gcf,'position',[440   159   802   639]);
hold on;
% for frame1 = 1:params.nSteps%size(movingBar.sceneRGB,3)
for frame1 = [1:size(movieStitch,3)-1]
    subplot(221);
    imagesc(squeeze(movingBar.sceneRGB(:,:,frame1,:)));
    caxis([0 1]);
    colormap gray; 
    
    subplot(222);
    for xPos = 3:numberElectrodesX-2
        for yPos = 1:numberElectrodesY
            hold on;
            fill(xh+electrodeArray.center(xPos,numberElectrodesY+1-yPos,1),yh+electrodeArray.center(xPos,numberElectrodesY+1-yPos,2),electrodeArray.activation(xPos,yPos,frame1))
        end
    end
    caxis([0 1]);
    
    
    subplot(223);    
%     imagesc((stimulusReconstructionHealthy(1:paramsRecHealthy.maxx,1:paramsRecHealthy.maxy,frame1)));
%     imagesc((stimulusReconstructionHealthy(1:sizeScene(2),1:sizeScene(2),frame1)));
%         imagesc(movrecons_on(:,:,frame1));

        imagesc(movieStitchHealthy(:,:,frame1-0));
%         imagesc(movieStitchHealthyDropout(:,:,frame1-20)');
     colormap gray
%      caxis([1*paramsRecHealthy.minR 1*paramsRecHealthy.maxR]);
%     caxis([1*paramsRecHealthy.minR 1*paramsRecHealthy.maxR]);
%     caxis([0 .5*paramsRecHealthy.maxR]);
%     caxis([0 .75*max(movrecons_on(:))]);
    caxis(caxishealthy);
    title('Healthy');
    
    subplot(224);    
%     imagesc((stimulusReconstruction(1:paramsRec.maxx,1:paramsRec.maxy,frame1)));
%     imagesc((stimulusReconstruction(1:sizeScene(2),1:sizeScene(2),frame1))');
% imagesc((stimulusReconstruction(1:200,1:200,frame1))');
% 
% % moviemat = movrecons_on;
%      colormap gray
% % %     caxis([.5*paramsRec.minR .5*paramsRec.maxR]);
%     caxis([1*paramsRec.minR 1*paramsRec.maxR]);

    imagesc(movieStitch(:,:,frame1-0));  
%         imagesc(movieStitchDropout(:,:,frame1-20-4)');

     colormap gray
%      caxis([1*paramsRecHealthy.minR 1*paramsRecHealthy.maxR]);
%     caxis([1*paramsRecHealthy.minR 1*paramsRecHealthy.maxR]);
%     caxis([0 .5*paramsRecHealthy.maxR]);
%     caxis([0 .75*max(movrecons_on(:))]);
    caxis(caxispros);
%     title('Healthy');
    title('Prosthetic - 0% dead');
%     pause(0.1);
drawnow

    F = getframe(h1);
    writeVideo(vObj,F);
end
end


close(vObj)
% 
% 
% % close all;
% % clear stimulusReconstruction stimulusReconstructionHealthy
% % name_str = ['gratingH_20Hz_width_' num2str(params.barWidth) '_freq_' num2str(freqL) '_onM_8_hz_ON_IMS1_' num2str(cputime*100) '.mat'];
% % save([path_str name_str])