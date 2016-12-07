% t_pixiumTestDirect
%
% Reconstructs a stimulus from direct electrode activation
%
% Stimulating the activations of retinal ganglion cells with an array of
% electrodes. The stimulus image is defined, the electrode array is
% generated, the electrode activations are computed, the RGC mosaics are
% generated, the RGC mosaic responses are computed and the stimulus is
% inferred from the RGC mosaic resposnes using a linear decoder.
%
% Outline of computation:
% 1. Load image/movie
% 2. Outer segment representation
% 3. Build electrode array
% 4. Compute electrode activations from image/movie
% 5. Build RGC array
% 6. Calculate RGC input - only one, change to input from multiple
% 7. Build RGC activation functions
% 8. Compute RGC activations/spikes
% 9. Invert representation to form image/movie
% 10. Tile over big FOV
%
% 12/2016 JRG (c) isetbio
%
%  3bd0152 % nov 19\

%% Initialize
clear;
% ieInit;
tic

%% Parameters to alter
clear electrodeArray

% Electrode size
% Set the size of implant pixels
electrodeArray.width = 15e-6; % meters
% electrodeArray.width = 140e-6; % meters

% % % % % % KNOBS
pulseFreq = 25;           % Hz, electrode pulse frequency
pulseDutyCycle = 1;       % Fraction of cycle pulse is on
irradianceFraction = .5;  % Fraction of maximum irradiance 

% Retinal patch eccentricity
patchEccentricity = 4; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.6*2/3;

% Stimulus length
nSteps = 520;

% %%% Grating subunit stimulus
%
% % params.barWidth = bwL;
% % iStim = ieStimulusGratingSubunit(params);
%
% % iStim = iStimC;
% absorptions = iStim.absorptions;
% % movingBar = iStim;
%
% nSteps = size(iStim.sceneRGB,3);
% params.nSteps = nSteps;
%
% size(iStim.sceneRGB)

percentDead = 0;

params.nSteps = nSteps;
params.row = 96;
params.col = 96;
params.fov = 1.6;

%% Loop for different resizing factors
% The hall movies captures a large field of view. In order to reconstruct
% this large-FOV stimulus with a small RGC mosaic, we must tile the mosaic
% over the stimulus a number of times. The mosaic is necessarily small
% because the training algorithm takes a long time.

rsFactor = 1;%[1 2 3 5 6]

%% Resize the hallway movie stimulus for tiling
rsFactor
tic

% Stimulus parameters
paramsStim.nsteps = 1;%nFrames;%size(testmovieshort,3);
paramsStim.timeInterval = 1/125;%0.001; % sec
paramsStim.expTime = 1/125;%0.001; % sec
nFrames = nSteps;

paramsStim.nsteps = nFrames;
paramsStim.fov = 8;
paramsStim.radius = 36/3*1e-6;
paramsStim.theta = 330;
paramsStim.side = 'left';
paramsStim.fov = fov;
blockctr = 0;

tic

%% Loop over each tiled mosaic
iblock = 1;%:rsFactor
jblock = 1;%:rsFactor

% Get iStim structure with os object with desired properties
blockctr = blockctr+1;
paramsStim.nsteps = 1;
testmovieshort = rand(96,96,2);
iStim = ieStimulusMovie(testmovieshort,paramsStim);

iStim.absorptions = iStim.sensor;
clear movingBar
movingBar = iStim;

%% Outer segment calculation
% There is no simulated outer segment, this identity outer segment acts as
% a pass-through holder of the stimulus intensity information.

% Input = RGB
os = osCreate('displayrgb');

% Get retinal patch properties for electrode array
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

movingBar.sceneRGB = testmovieshort((iblock-1)*96+[1:96],(jblock-1)*96+[1:96],:)
% movingBar.sceneRGB = (contrastElectrode)*(movingBar.sceneRGB - 0.5)+0.5;
os = osSet(os, 'rgbData', movingBar.sceneRGB);

sceneRGB_Healthy = (irradianceFraction)*(movingBar.sceneRGB - 0.5)+0.5;
osHealthy = os;
osHealthy = osSet(osHealthy, 'rgbData', sceneRGB_Healthy);

retinalPatchSize = osGet(os,'size');

% Electrode array properties
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

th = (0:1/6:1)'*2*pi;
xh = electrodeArray.width/2*cos(th);
yh = electrodeArray.width/2*sin(th);

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
        
        % sizeES = size(electrodeStimulus);
        % electrodeArray.activation(xPos,yPos,frame) = min([ mean(electrodeStimulus(:,1:floor(sizeES(2)/2))) mean(electrodeStimulus(:,ceil(sizeES(2)/2):sizeES(2)))]);
    end
end

electrodeArray.activation(:) = 0;

for telec = 1:numberElectrodesX
    electrodeArray.activation(telec,:,5*(telec-1)+1:5*telec) = 1;
end
% for telec = 1:4*numberElectrodesX
%     electrodeArray.activation(telec,:,5*(floor((telec+3)/4)-1)+1:5*floor((telec+1)/4)) = 1;
% end
%% Add X Hz spiking of stimulus

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

%% Build RGC array for healthy retina
clear paramsIR innerRetinaHealthy
load('C:\Users\James\Documents\MATLAB\github\RGC-Reconstruction\dat\pixium\mosaicAll_pix_ns.mat')

%% Calculate RGC input
% Weight electrode activation by Gaussian as a function of distance between
% centers
offFlag = 0; nTileRows = 1; nTileCols = 1; mosaicOffset = zeros(4,1);
innerRetinaInput = irActivationFromElectrode(innerRetina, electrodeArray, retinalPatchWidth, metersPerPixel, nTileRows, nTileCols, mosaicOffset, params, offFlag);

%% Build RGC activation functions
innerRetinaThreshold = irGetThreshold(innerRetina);

%% Compute RGC activations
innerRetina =  irGetLinearRespElectrode(innerRetina, 10*innerRetinaInput, innerRetinaThreshold, params);

%% Compute RGC spiking
numberTrials = 1;
for tr = 1:numberTrials
    innerRetina = irComputeSpikes(innerRetina,'coupling',false);
end
clear innerRetina1
innerRetina1 = innerRetina;
toc
%% Do optimal reconstruction

pOpt.innerRetina = innerRetina;
pOpt.percentDead = 0;
pOpt.numbins = 4;
pOpt.filterFile = 'pixium\pix1_long_filter_nsBig_100hz_4st';

[movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalReconSingle(pOpt);

movieComb = 255*irradianceFraction*pulseDutyCycle*ieScale(movrecons_on_offHealthy);
szFrames = size(electrodeArray.activation,3);
eAact = zeros(rsFactor*96,rsFactor*96,szFrames);
for ii = 1:szFrames
    eAact(:,:,ii) = imresize(electrodeArray.activation(:,:,ii),[rsFactor*96,rsFactor*96])';
end
movieComb(:,rsFactor*96+[1:rsFactor*96],1:size(eAact,3)) = 255*irradianceFraction*pulseDutyCycle*ieScale(eAact(:,:,1:size(eAact,3)));% testmovieshort(:,:,shiftval+1:567+1);
movieComb(:,:,end) = 255;
% figure; ieMovie(movieComb);
%% Save for tiling
% save([reconstructionRootPath '\dat\pixium\dir_dec5_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetina','movrecons_on_offHealthy');
toc

% %% Save as output with original movie
fig=figure;
% set(fig,'position',[    624         437        1018         541]);
set(fig,'position',[624   704   537   274]);
aviobj = avifile([reconstructionRootPath '\dat\pixium\prosthesis_recon_' num2str(rsFactor) '_directBar.avi'])
aviobj.Fps = 30;

for k=1:size(movieComb,3)-0
    % imagesc(movieRecon(:,:,k)); colormap gray;
    
    image(uint8(movieComb(:,:,k))); colormap gray; axis image
    caxis([0 255]);
    F = getframe(fig);
    aviobj = addframe(aviobj,F);
end
close(fig)
aviobj = close(aviobj);

%%
clear movieComb movieRecon testmovieshort vidFrame
toc