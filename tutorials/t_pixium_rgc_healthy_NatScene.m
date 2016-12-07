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
% Add 5 Hz stimulation pulsesx

%  3bd0152
% 3bd015266549daddbb470052bcecfe5433154838 % nov 19


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
electrodeArray.width = 2*15e-6; % meters
% electrodeArray.width = 140e-6; % meters
 
% Retinal patch eccentricity
patchEccentricity = 4; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.6*2/3;

% Stimulus length
nSteps = 520;

% Activation curve

% Spatial activation of electrode

% Electrode PWM 

percentDead = 0;
%% Load image

% clear params
% % One frame of a moving bar stimulus
% % Set parameters for size
% % params.nSteps = nSteps;
% params.barWidth = 20;
% params.row = 96;
% params.col = 96;
% params.fov = fov;
% % params.freq = freqL; % Hz grating frequency
% % % params.vfov = 0.7;
% params.startFrames = 0;
% params.endFrames = 0;
% movingBar = ieStimulusBarSensor(params);
% iStim = movingBar;
% tuningWoffElec = 1;
% tuningWoffHealthy = 1;
% 
pulseFreq = 25; % Hz, electrode pulse frequency
% 
contrastHealthy = 1;
contrastElectrode = 1;
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


%%


% Length of WN movie is 1200, take nFrames to limit natural movie to same length
% nFrames = 3600; 
% testmovieshort = double(testmovie.matrix(:,:,1:nFrames)); 
% load('EJLPhosphene/dat/stimuli/hallMovie.mat')
% testmovieshort = vidFrame; clear vidFrame;

% load('C:\Users\James\Documents\MATLAB\github\EJLPhosphene\dat\stimuli\hallMovie288.mat','hallMovieResize');

for rsFactor = 4%[1 2 3 5 6]
rsFactor
tic
load('C:\Users\James\Documents\MATLAB\github\EJLPhosphene\dat\stimuli\hallMovie.mat')
% load([reconstructionRootPath '\dat\stimuli\hallMovie.mat'])
szFrames = size(vidFrame,3);
hallMovieResize = zeros(rsFactor*96,rsFactor*96,szFrames);
for ii = 1:szFrames
    hallMovieResize(:,:,ii) = imresize(vidFrame(:,:,ii),[rsFactor*96,rsFactor*96]);
end

testmovieshort = (255*ieScale(hallMovieResize)); clear hallMovieResize;
% testmovieshort = uint8(255*ieScale(hallMovieResize)); clear hallMovieResize;
% Generate display, scene, oi, sensor
paramsStim.nsteps = 1;%nFrames;%size(testmovieshort,3);
%  Bipolar filter is setfor 0.001 sec, so it needs to be 0.001
paramsStim.timeInterval = 1/125;%0.001; % sec
paramsStim.expTime = 1/125;%0.001; % sec
nFrames = nSteps;
% % For 2013-08-19-6
% r = 12 mm = 36 degs;
% theta = 330 degs;
% side = 'left';
paramsStim.nsteps = nFrames;
paramsStim.fov = 8;
paramsStim.radius = 36/3*1e-6;
paramsStim.theta = 330;
paramsStim.side = 'left';
paramsStim.fov = fov;
% iStim = ieStimulusMovie(testmovieshort(:,:,1:nFrames),paramsStim);
blockctr = 0;
tic
for iblock = 1:rsFactor
    for jblock = 1:rsFactor
       blockctr = blockctr+1;
% iStim = ieStimulusMovie(testmovieshort(96+[1:96],96+[1:96],1:nFrames),paramsStim);
paramsStim.nsteps = 10;
iStim = ieStimulusMovie(testmovieshort((iblock-1)*96+[1:96],(jblock-1)*96+[1:96],1:10),paramsStim);

% paramsStim.nsteps = 10;
% iStim = ieStimulusMovieCMosaic(testmovieshort(96+[1:96],96+[1:96],1:10),paramsStim);

iStim.absorptions = iStim.sensor;
clear movingBar
movingBar = iStim;
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

movingBar.sceneRGB = testmovieshort((iblock-1)*96+[1:96],(jblock-1)*96+[1:96],:)
% movingBar.sceneRGB = (contrastElectrode)*(movingBar.sceneRGB - 0.5)+0.5;
os = osSet(os, 'rgbData', movingBar.sceneRGB);

sceneRGB_Healthy = (contrastHealthy)*(movingBar.sceneRGB - 0.5)+0.5;
osHealthy = os;
osHealthy = osSet(osHealthy, 'rgbData', sceneRGB_Healthy);


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


% filenameRGC = ['/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/WNstim_response_OffMidget_RGC.mat'];
filenameRGC = ['C:\Users\James\Documents\GitHub\offMidget1\WNstim_response_OffMidget_RGC.mat'];
load(filenameRGC);
innerRetinaHealthy3 = innerRetina;
% for loadind = 2:4
%     clear innerRetina
%     filenameRGC = ['/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/WNstim_response_OffMidget_RGC.mat'];
%     load(filenameRGC);
%     innerRetinaHealthy3.mosaic{loadind} = innerRetina.mosaic{1};
%     
% end

% filenameRGC = ['/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/WNstim_response_OnMidget_RGC.mat'];
filenameRGC = ['C:\Users\James\Documents\GitHub\onMidget1\WNstim_response_OnMidget_RGC.mat'];
load(filenameRGC);
innerRetinaHealthy4 = innerRetina;
% for loadind = 2:4
%      clear innerRetina
%     filenameRGC = ['/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/WNstim_response_OnMidget_RGC.mat'];
%     load(filenameRGC);
%     innerRetinaHealthy4.mosaic{loadind} = innerRetina.mosaic{1};
%     
% end

% rdt = RdtClient('isetbio');
% rdt.crp('/resources/data/rgc/pixium/');
% data = rdt.readArtifact('WNstim_response_OffParasol_RGC', 'type', 'mat');
filenameRGC = ['C:\Users\James\Documents\GitHub\May26_offBig2\WNstim_response_OffParasol_RGC.mat'];
load(filenameRGC);
innerRetinaHealthy2 = innerRetina; % clear innerRetina;
% innerRetinaHealthy2 = data.innerRetina; % clear innerRetina;
% for loadind = 2:4
%     clear data
%     data = rdt.readArtifact('WNstim_response_OffParasol_RGC', 'type', 'mat');
%     innerRetinaHealthy2.mosaic{loadind} = data.innerRetina.mosaic{1};
%     
% end

filenameRGC = ['C:\Users\James\Documents\GitHub\may26_onBig2\WNstim_response_OnParasol_RGC.mat'];
load(filenameRGC);
% data = rdt.readArtifact('WNstim_response_OnParasol_RGC', 'type', 'mat');
% innerRetinaHealthy = data.innerRetina;
innerRetinaHealthy = innerRetina;
% for loadind = 2:4
%     clear data
%     data = rdt.readArtifact('WNstim_response_OnParasol_RGC', 'type', 'mat');
%     innerRetinaHealthy.mosaic{loadind} = data.innerRetina.mosaic{1};
%     
% end

innerRetinaHealthy = irCompute(innerRetinaHealthy,osHealthy,'coupling',false);
innerRetinaHealthy2 = irCompute(innerRetinaHealthy2,osHealthy,'coupling',false);
innerRetinaHealthy3 = irCompute(innerRetinaHealthy3,osHealthy,'coupling',false);
innerRetinaHealthy4 = irCompute(innerRetinaHealthy4,osHealthy,'coupling',false);

%%

% clear movrecons_on_off movieStitch
[movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalRecon(innerRetinaHealthy, innerRetinaHealthy2, innerRetinaHealthy3, innerRetinaHealthy4, percentDead);


%%

save(['C:\Users\James\Documents\MATLAB\github\EJLPhosphene\dat\tile\ws_nov24_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetinaHealthy', 'innerRetinaHealthy2', 'innerRetinaHealthy3', 'innerRetinaHealthy4','movrecons_on_offHealthy');
toc
    end
end

%%
movieRecon = zeros(rsFactor*96,rsFactor*96,591);
blockctr = 0;
tic
for iblock = 1:rsFactor
    for jblock = 1:rsFactor
       blockctr = blockctr+1;
       
%        load(['ws_nov24_block' num2str(blockctr) '.mat']);%, 'innerRetinaHealthy', 'innerRetinaHealthy2', 'innerRetinaHealthy3', 'innerRetinaHealthy4');
load(['C:\Users\James\Documents\MATLAB\github\EJLPhosphene\dat\tile\ws_nov24_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat']);

       movieRecon((iblock-1)*96+[1:96],(jblock-1)*96+[1:96],:) = movrecons_on_offHealthy{1}(:,:,1:591);
       clear movrecons_on_offHealthy
    end
end

save(['C:\Users\James\Documents\MATLAB\github\EJLPhosphene\dat\nov29_movRecon_' num2str(rsFactor) '.mat'],'movieRecon');

% figure; ieMovie(movieRecon);

%%
fig=figure;
set(fig,'position',[    624         437        1018         541]);
aviobj = avifile(['C:\Users\James\Documents\MATLAB\github\EJLPhosphene\dat\prosthesis_recon_' num2str(rsFactor) '_30fps_both2.avi'])
aviobj.Fps = 30;
shiftval = 13;
movieComb = 355*ieScale(movieRecon(:,:,1:591-shiftval+1));
movieComb(:,rsFactor*96+[1:rsFactor*96],:) = testmovieshort(:,:,shiftval+1:591+1);
    
for k=1:size(movieRecon,3)-60
    % imagesc(movieRecon(:,:,k)); colormap gray;

    imagesc(movieComb(:,:,k)); colormap gray; axis image
    F = getframe(fig);
    aviobj = addframe(aviobj,F);
end
close(fig)
aviobj = close(aviobj);
 
%%
clear movieComb movieRecon testmovieshort vidFrame
toc
end