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
saveFile = 'pix1';

% bwL = 20%[ 32 50]
% freqL = 5%[5 8 2 12]

pulseFreq = 100; % Hz, electrode pulse frequency
%% Parameters to alter
clear electrodeArray
% Electrode size
% Set the size of implant pixels
electrodeArray.width = 1*25e-6; % meters
% electrodeArray.width = 140e-6; % meters
 
% Retinal patch eccentricity
patchEccentricity = 4; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.6*2/3;

% Activation curve

% Spatial activation of electrode

% Electrode PWM 

percentDead = 0;
%% Load image

%% Parameters to alter

% Retinal patch eccentricity
patchEccentricity = 12; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 2*1.6;

% Stimulus length = nSteps*nBlocks;
nSteps = 1;
nBlocks = 80;

%% Load image
clear params
% One frame of a moving bar stimulus
% Set parameters for size
params.nSteps = nSteps;
params.row = 200;
params.col = 200;
params.fov = fov;
% % params.vfov = 0.7;

%%


iStim = ieStimulusBinaryWhiteNoise(params);
absorptions = iStim.sensor;
whiteNoise = iStim;

nSteps = 1;

% Input = RGB
os = osCreate('displayrgb');

sceneSize = sceneGet(whiteNoise.scene,'size');
retinalPatchWidth = sensorGet(whiteNoise.sensor,'width','m');
% retinalPatchHeight = sensorGet(whiteNoise.absorptions,'height','m');
retinalPatchHeight = (sceneSize(1)/sceneSize(2))*retinalPatchWidth;

% % % coneSpacing = scene.wAngular*300
% coneSpacing = sensorGet(sensor,'dimension','um');
os = osSet(os, 'patchSize', retinalPatchWidth);

% timeStep = sensorGet(whiteNoise.sensor,'time interval','sec');
timeStep = (1/125)/1;
os = osSet(os, 'timeStep', timeStep);

% os = osSet(os, 'rgbData', whiteNoise.sceneRGB);
% os = osCompute(absorptions);


% load([ reconstructionRootPath  '\dat\movsm_' num2str(blockNum) '.mat'],'movsm');
%      load([ reconstructionRootPath  '/dat/imagenetBlocks/movsm_' num2str(blockNum) '.mat'],'movsm');
%           load([ reconstructionRootPath  '\dat\imagenetBlocks\movsm_' num2str(blockNum) '.mat'],'movsm');
% natScenes = movsm(1:96,1:96,randperm(nSteps));
os = osSet(os, 'rgbData', rand(200,200,2));

%% Build RGC array

% if isempty(mosaicFile)
    % clear paramsIR innerRetina
    paramsIR.name    = 'Macaque inner retina 1'; % This instance
    paramsIR.eyeSide   = 'left';   % Which eye
    paramsIR.eyeRadius = 2.5;        % Radius in mm
    paramsIR.eyeAngle  = 90;       % Polar angle in degrees
    
    model   = 'LNP';    % Computational model
%     innerRetina = irCreate(os,paramsIR);
%     innerRetina = rgcMosaicCreate(innerRetina,'type','onParasol','model',model);
%     innerRetina = rgcMosaicCreate(innerRetina,'type','offParasol','model',model);
%     innerRetina = rgcMosaicCreate(innerRetina,'type','offMidget','model',model);
%     innerRetina = rgcMosaicCreate(innerRetina,'type','onMidget','model',model);
    
    % innerRetina = rgcMosaicCreate(innerRetina,'type','sbc','model',model);
    
    % innerRetina.mosaic{1}.mosaicSet('numberTrials',1);
    % innerRetina.mosaic{2}.mosaicSet('numberTrials',1);
    % innerRetina.mosaic{3}.mosaicSet('numberTrials',1);
    % innerRetina.mosaic{4}.mosaicSet('numberTrials',1);
    % innerRetina.mosaic{5}.mosaicSet('numberTrials',1);
    
    % irPlot(innerRetina,'mosaic');
    
%     mosaicFile = ['mosaicAll_' num2str(round(cputime*100))];
%     filenameRGC = [phospheneRootPath '/dat/' mosaicFile '.mat'];
%     save(filenameRGC, 'innerRetina');
    
    mosaicFile = ['mosaicAll_8372855'];
    filenameRGC = [phospheneRootPath '/dat/' mosaicFile '.mat'];
    load(filenameRGC);

% else
%     
%     % filenameRGC = [reconstructionRootPath '\dat\mosaic_all.mat'];
%     % filenameRGC = [reconstructionRootPath '/dat/mosaic_all_bertha_ns0.mat'];
%     filenameRGC = [reconstructionRootPath '/dat/' mosaicFile '.mat'];
% end


for blockNum =3:14

%%% Grating subunit stimulus
tic
iStim = ieStimulusBinaryWhiteNoise(params);
absorptions = iStim.sensor;
whiteNoise = iStim;
%% Show raw stimulus for osIdentity
% % figure;
% % for frame1 = 1:size(whiteNoise.sceneRGB,3)
% %     imagesc(squeeze(whiteNoise.sceneRGB(:,:,frame1,:)));
% %     colormap gray; drawnow;
% % end
% % % close;

%% Outer segment calculation
% There is no simulated outer segment, this identity outer segment acts as
% a pass-through holder of the stimulus intensity information.

nSteps = 12000;

% Input = RGB
os = osCreate('displayrgb');

sceneSize = sceneGet(whiteNoise.scene,'size');
retinalPatchWidth = sensorGet(whiteNoise.sensor,'width','m');
% retinalPatchHeight = sensorGet(whiteNoise.absorptions,'height','m');
retinalPatchHeight = (sceneSize(1)/sceneSize(2))*retinalPatchWidth;

% % % coneSpacing = scene.wAngular*300
% coneSpacing = sensorGet(sensor,'dimension','um');
os = osSet(os, 'patchSize', retinalPatchWidth);

% timeStep = sensorGet(whiteNoise.sensor,'time interval','sec');
timeStep = (1/125)/1;
os = osSet(os, 'timeStep', timeStep);

% os = osSet(os, 'rgbData', whiteNoise.sceneRGB);
% os = osCompute(absorptions);


% load([ reconstructionRootPath  '\dat\movsm_' num2str(blockNum) '.mat'],'movsm');
     load([ reconstructionRootPath  '/dat/imagenetBlocks200/movsm_' num2str(blockNum) '.mat'],'movsm');
%           load([ reconstructionRootPath  '\dat\imagenetBlocks\movsm_' num2str(blockNum) '.mat'],'movsm');
natScenes = movsm(:,:,randperm(nSteps));
os = osSet(os, 'rgbData', double(natScenes));

% % Plot the photocurrent for a pixel
% osPlot(os,absorptions);

retinalPatchSize = osGet(os,'size');
% os = osCompute(absorptions);

% % Plot the photocurrent for a pixel
% osPlot(os,absorptions);

metersPerPixel = retinalPatchWidth/retinalPatchSize(2);

retinalPatchSize = osGet(os,'size');
numberElectrodesX = floor(retinalPatchWidth/electrodeArray.width)+4;
numberElectrodesY = floor(retinalPatchWidth/electrodeArray.width)+4;
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
        electrodeStimulus = squeeze(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2,:,:));
        electrodeArray.activation(xPos,yPos,:) = mean(RGB2XWFormat(electrodeStimulus));
        
        %             sizeES = size(electrodeStimulus);
        %             electrodeArray.activation(xPos,yPos,frame) = min([ mean(electrodeStimulus(:,1:floor(sizeES(2)/2))) mean(electrodeStimulus(:,ceil(sizeES(2)/2):sizeES(2)))]);
        
        % imagesc(electrodeStimulus); title(sprintf('%2.2f',mean(electrodeStimulus(:))));
    end
end


% Show plot

% Apply Gaussian


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

%%
load(filenameRGC);

%% Calculate RGC input
% Weight electrode activation by Gaussian as a function of distance between
% centers
% params.nSteps = paramsStim.nsteps-8;
% On midget
offFlag = [0 1 1 0]; nTileRows = 1; nTileCols = 1; mosaicOffset = zeros(4,1);
innerRetinaInput = irActivationFromElectrode(innerRetina, electrodeArray, retinalPatchWidth, metersPerPixel, nTileRows, nTileCols, mosaicOffset, params, offFlag);

%% Build RGC activation functions

innerRetinaThreshold = irGetThreshold(innerRetina);

%% Compute RGC activations

innerRetina =  irGetLinearRespElectrode(innerRetina, 1*innerRetinaInput, innerRetinaThreshold, params);
%% Compute RGC spiking
numberTrials = 1;
for tr = 1:numberTrials
    innerRetina = irComputeSpikes(innerRetina,'coupling',false);
end
% clear innerRetina1
% innerRetina1 = innerRetina;
toc
%%
tic   
    spikesout  = RGB2XWFormat(mosaicGet(innerRetina.mosaic{1},'spikes'));
    spikesout2 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{2},'spikes'));
    spikesout3 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{3},'spikes'));
    spikesout4 = RGB2XWFormat(mosaicGet(innerRetina.mosaic{4},'spikes'));
    
    timeBins = max([size(spikesout,2) size(spikesout2,2) size(spikesout3,2) size(spikesout4,2)]);
    
    spikesoutsm = zeros(size(spikesout,1)+ size(spikesout2,1)+size(spikesout3,1)+size(spikesout4,1), timeBins,'uint8');
    spikesoutsm(1:size(spikesout,1) ,1:size(spikesout,2) ) = spikesout;
    spikesoutsm(size(spikesout,1)+[1:size(spikesout2,1)],1:size(spikesout2,2) ) = spikesout2;
    
    spikesoutsm(size(spikesout,1)+size(spikesout2,1)+[1:size(spikesout3,1)] ,1:size(spikesout3,2) ) = spikesout3;
    spikesoutsm(size(spikesout,1)+size(spikesout2,1)+size(spikesout3,1)+[1:size(spikesout4,1)] ,1:size(spikesout4,2) ) = spikesout4;
  
%     whiteNoiseSmall = uint8(squeeze(whiteNoise.sceneRGB(:,:,:,1)));
%     whiteNoiseSmall = iStim.sceneRGB;
    whiteNoiseSmall = natScenes;
    
    % spikesoutsm = uint8(spikesoutmat);
    
    % filename1 = [reconstructionRootPath '\dat\NSstim_response_overlap0_block_' num2str(blockNum) '.mat'];
    % filename1 = [reconstructionRootPath '/dat/nsResponses/NSstim_response_betha_ns0_block_' num2str(blockNum) '.mat'];
%     filename1 = [reconstructionRootPath '/dat/nsResponses/' saveFile '_block_' num2str(blockNum) '.mat'];
    filename1 = [phospheneRootPath '/dat/pixium25/' saveFile '_nsBig_100Hz_block_' num2str(blockNum) '_' mosaicFile '.mat'];
    
    save(filename1, 'spikesoutsm','whiteNoiseSmall');
    toc
    %%
    
    clear spikesout spikesout2 spikesout3 spikesout4 spikesoutmat whiteNoise whiteNoiseSmall spikesoutsm fullStimulus % iStim
    clear innerRetinaInput innerRetinaInput2 innerRetinaInput3 innerRetinaInput4 innerRetina
    close all;
end
toc
% figure; imagesc(reshape(sum(abs(filterMat),1),[96 96]));
%%
clear spikesoutSum
for blockSp = 1:1000-1
    spikesoutSum(:,blockSp) = sum(spikesoutsm(:,10*(blockSp-1)+1:10*blockSp),2);
end
spikesoutSum(:,1000) = zeros(size(spikesoutsm,1),1);
% spikesoutSum = spikesoutsm;
rfout = single(RGB2XWFormat(whiteNoiseSmall))*single(spikesoutSum');
figure; 
cistart = 100;
for ci = 1:25
    subplot(5,5,ci);
    imagesc(reshape(rfout(:,ci+cistart),[200 200]))
end
% 
% %%
% 
% % figure; plot(squeeze(innerRetina2.mosaic{1}.responseLinear(1,1,:)))
% 
% sp1 = innerRetina3.mosaic{1}.responseSpikes{15,15};
% 
% strf = zeros(96,96,15);
% for sp = 1:length(sp1)-1
%     fr = round(innerRetina3.mosaic{1}.responseSpikes{15,15}(sp));
%     strf = strf + iStim.sceneRGB(:,:,fr:fr+14) - 0.5;
% %     strf = strf + whiteNoiseSmall(
%     
% end
%     
% 
% figure; imagesc(strf(:,:,2))
%     
% figure; 
% for fr = 1:15
%     subplot(4,4,fr)
%     imagesc(strf(:,:,fr))
% end
% 
