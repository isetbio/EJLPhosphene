% t_pixiumTestFullField
% 
% Reconstructs a full-field stimulus from direct electrode activation.
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

% Retinal patch eccentricity
patchEccentricity = 4; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.6*2/3;

% % % % % % KNOBS
pulseFreq = 25;           % Hz, electrode pulse frequency
pulseDutyCycle = .5;       % Fraction of cycle pulse is on
irradianceFraction = .5;  % Fraction of maximum irradiance 

% Stimulus length
nSteps = 520;

percentDead = 0;

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

params.nSteps = nSteps;
params.row = 96;
params.col = 96;
params.fov = 1.6;

%% Loop for different resizing factors
% The hall movies captures a large field of view. In order to reconstruct
% this large-FOV stimulus with a small RGC mosaic, we must tile the mosaic
% over the stimulus a number of times. The mosaic is necessarily small
% because the training algorithm takes a long time.

rsFactor = 1%[1 2 3 5 6]

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

sceneRGB_Healthy = (1)*(movingBar.sceneRGB - 0.5)+0.5;
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

%% Full field direct stimulation

electrodeArray.activation(:) = 0;

% % moving bar
% for telec = 1:25
%     electrodeArray.activation(telec,:,5*(telec-1)+1:5*telec) = 1;
% end

% full field
intervalSteps = 24;
telec = 2;
electrodeArray.activation(:,:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

% half field 1
telec = telec+2;
electrodeArray.activation(:,1:round(numberElectrodesY/2)-1,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

% % half field 2
telec = telec+2;
electrodeArray.activation(:,round(numberElectrodesY/2):end,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

% % half field 3
% telec = telec+2;
% electrodeArray.activation(1:12,:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

% % half field 4
% telec = telec+2;
% electrodeArray.activation(13:end,:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

%  quarter field 1
telec = telec+2;
electrodeArray.activation(round(numberElectrodesX/2):end,1:round(numberElectrodesY/2)-1,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

%  quarter field 2
telec = telec+2;
electrodeArray.activation(1:round(numberElectrodesX/2)-1,round(numberElectrodesY/2):end,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

% subsets of six at four points
% point 1
telec = telec+2;
tint = [intervalSteps*(telec-1)+1:intervalSteps*telec ...
        intervalSteps*(telec+1):intervalSteps*(telec+2) ...
        intervalSteps*(telec+3):intervalSteps*(telec+4)];
    
centerv = [round(numberElectrodesX/3) round(numberElectrodesY/3)]; hexcenter = [];
cctr = 0; for c1 = -1:1; for c2 = -1:1; cctr = cctr+1; hexcenter(cctr,:) = [centerv(1)+c1 centerv(2)+c2]; end; end;
electrodeArray.activation(hexcenter(:,1),hexcenter(:,2),tint) = 1;
% electrodeArray.activation([1:2 end-1:end],:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% electrodeArray.activation(:,[1:2 end-1:end],intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% point 2
centerv = [round(numberElectrodesX/3) round(2*numberElectrodesY/3)]; hexcenter = [];
cctr = 0; for c1 = -1:1; for c2 = -1:1; cctr = cctr+1; hexcenter(cctr,:) = [centerv(1)+c1 centerv(2)+c2]; end; end;
electrodeArray.activation(hexcenter(:,1),hexcenter(:,2),tint) = 1;
% electrodeArray.activation([1:2 end-1:end],:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% electrodeArray.activation(:,[1:2 end-1:end],intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% point 3
centerv = [round(2*numberElectrodesX/3)-1 round(numberElectrodesY/3)]; hexcenter = [];
cctr = 0; for c1 = -1:1; for c2 = -1:1; cctr = cctr+1; hexcenter(cctr,:) = [centerv(1)+c1 centerv(2)+c2]; end; end;
electrodeArray.activation(hexcenter(:,1),hexcenter(:,2),tint) = 1;
% electrodeArray.activation([1:2 end-1:end],:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% electrodeArray.activation(:,[1:2 end-1:end],intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% point 4
centerv = [round(2*numberElectrodesX/3)-1 round(2*numberElectrodesY/3)]; hexcenter = [];
cctr = 0; for c1 = -1:1; for c2 = -1:1; cctr = cctr+1; hexcenter(cctr,:) = [centerv(1)+c1 centerv(2)+c2]; end; end;
electrodeArray.activation(hexcenter(:,1),hexcenter(:,2),tint) = 1;
% electrodeArray.activation([1:2 end-1:end],:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% electrodeArray.activation(:,[1:2 end-1:end],intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;


%%%% Border activation
xborder = 2; yborder = 1;
electrodeArray.activation([1:xborder end-xborder+1:end],:,:) = .25;
electrodeArray.activation(:,[1:yborder end-yborder+1:end],:) = .25;

% telec = telec+2;
electrodeArray.activation(:,:,intervalSteps*(telec+4):intervalSteps*(telec+5)) = 0;

% figure; ieMovie(electrodeArray.activation,'FrameRate',10);
%% Add X Hz spiking of stimulus

% Right now the electrode sampling is at 0.01 s = 100 Hz
% Downsample to get 5 Hz
szAct = size(electrodeArray.activation);
electrodeArray.activationDS = zeros(szAct);
for iSample = 1:szAct(3)
    if mod(iSample,100/pulseFreq)==0
        electrodeArray.activationDS(:,:,iSample) = electrodeArray.activation(:,:,iSample);
%         electrodeArray.activationDSoff(:,:,iSample) = 1-electrodeArray.activation(:,:,iSample);
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
if isunix || ismac
    load([phospheneRootPath '/dat/mosaicAll_pix_ns.mat'])
else    
    load([phospheneRootPath '\dat\mosaicAll_pix_ns.mat'])
end

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
toc
%% Do optimal reconstruction

pOpt.innerRetina = innerRetina;
pOpt.percentDead = 0;
pOpt.numbins = 4;
pOpt.filterFile = 'pix1_long_filter_nsBig_100hz_4st';

[movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalReconSingle(pOpt);

% figure; ieMovie(movrecons_on_offHealthy);

shiftval = 13;
movieComb = 255*irradianceFraction*pulseDutyCycle*ieScale(movrecons_on_offHealthy);

szFrames = size(electrodeArray.activation,3);
eAact = zeros(rsFactor*96,rsFactor*96,szFrames);
for ii = 1:szFrames
    eAact(:,:,ii) = imresize(electrodeArray.activation(:,:,ii),[rsFactor*96,rsFactor*96])';
end

movieComb(:,rsFactor*96+[1:rsFactor*96],:) = 255*irradianceFraction*pulseDutyCycle*ieScale(eAact(:,:,1:size(movrecons_on_offHealthy,3)));% testmovieshort(:,:,shiftval+1:567+1);
% figure; ieMovie(movieComb);
%% Save for tiling
% save([reconstructionRootPath '\dat\pixium\dir_dec5_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetina','movrecons_on_offHealthy');
toc

% %% Save as output with original movie
fig=figure;
% set(fig,'position',[    624         437        1018         541]);
set(fig,'position',[624   704   537   274]);
if ismac || isunix    
    aviobj = avifile([phospheneRootPath '/dat/prosthesis_recon_' num2str(rsFactor) '_directFullField.avi'])
else
    aviobj = avifile([phospheneRootPath '\dat\prosthesis_recon_' num2str(rsFactor) '_directFullField.avi'])
end
aviobj.Fps = 30;

for k=1:size(movieComb,3)-0
    % imagesc(movieRecon(:,:,k)); colormap gray;
    
    image(movieComb(:,:,k)); colormap gray; axis image
    caxis([0 255]);
    F = getframe(fig);
    aviobj = addframe(aviobj,F);
end
close(fig)
aviobj = close(aviobj);

%%
clear movieComb movieRecon testmovieshort vidFrame
toc
