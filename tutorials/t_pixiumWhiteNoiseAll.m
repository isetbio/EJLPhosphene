% t_pixiumWhiteNoise
% 
% Run binary white noise through the RGC array set up for the Pixium
% prosethesis.


%% Initialize
clear;
% ieInit;

%% Parameters to alter

% Retinal patch eccentricity
patchEccentricity = 12; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.6;

% Stimulus length
nSteps = 2400;

% Activation curve

% Spatial activation of electrode

% Electrode PWM 


%% Load image
clear params
% One frame of a moving bar stimulus
% Set parameters for size
params.nSteps = nSteps;
params.row = 96;
params.col = 96;
params.fov = fov;
% % params.vfov = 0.7;

%%% Grating subunit stimulus

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

os = osSet(os, 'rgbData', whiteNoise.sceneRGB);
% os = osCompute(absorptions);

% % Plot the photocurrent for a pixel
% osPlot(os,absorptions);

retinalPatchSize = osGet(os,'size');

%% Build RGC array

clear paramsIR innerRetina
paramsIR.name    = 'Macaque inner retina 1'; % This instance
paramsIR.eyeSide   = 'left';   % Which eye
paramsIR.eyeRadius = 4;        % Radius in mm
paramsIR.eyeAngle  = 90;       % Polar angle in degrees

model   = 'LNP';    % Computational model
innerRetina = irCreate(os,paramsIR);
% innerRetina = rgcMosaicCreate(innerRetina,'type','offMidget','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','offMidget','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','onParasol','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','offParasol','model',model);
% innerRetina = rgcMosaicCreate(innerRetina,'type','sbc','model',model);
% innerRetina.mosaic{1}.mosaicSet('numberTrials',1);
% irPlot(innerRetina,'mosaic');
% % figure;

% innerRetina = irSet(innerRetina,'numberTrials',1);

% filenameRGC = ['/Users/james/Documents/MATLAB/isetbio misc/optimal linear decoder/May25_on2/WNstim_response_OffParasol_RGC_may16.mat'];

% filenameRGC = ['C:\Users\James\Documents\GitHub\offMidget1\WNstim_response_OffMidget_RGC.mat'];
% filenameRGC = ['C:\Users\James\Documents\GitHub\onMidget1\WNstim_response_OnMidget_RGC.mat'];

filenameRGC = ['C:\Users\James\Documents\GitHub\May26_onBig2\WNstim_response_OnParasol_RGC.mat'];

% filenameRGC = ['C:\Users\James\Documents\GitHub\May26_offBig1\WNstim_response_OffParasol_RGC.mat'];

% save(filenameRGC, 'innerRetina');


% load('/Users/james/documents/MATLAB/isetbio misc/optimal linear decoder/WNstim_response_OnParasol_RGC_may10.mat');

% 
load ../dat/movie_onMidget.mat

for blockNum = 1%:100%1:200

% clear psthNorm spikesout spikesoutM spikesoutsm whiteNoiseSmall whiteNoise iStim absorptions innerRetina

whiteNoiseSmall = reshape(stim(:,(blockNum-1)*1*2400+1:blockNum*2400),96,96,1*2400);
% clear stim
blockNum
%%% Grating subunit stimulus
load(filenameRGC, 'innerRetina');
% iStim = ieStimulusBinaryWhiteNoise(params);
% absorptions = iStim.sensor;
% whiteNoise = iStim;

% filename1 = ['C:\Users\James\Documents\GitHub\onMidget1\WNstim_response_OnMidget_block_' num2str(blockNum) '.mat'];
% filename1 = ['C:\Users\James\Documents\GitHub\offMidget1\WNstim_response_OffMidget_block_' num2str(blockNum) '.mat'];

% filename1 = ['C:\Users\James\Documents\GitHub\offParasol1_all\WNstim_response_OffParasol_block_' num2str(blockNum) '.mat'];
% load(filename1); clear spikesoutsm;
whiteNoise.sceneRGB = double(whiteNoiseSmall);

os = osSet(os, 'rgbData', whiteNoise.sceneRGB);

innerRetina = irCompute(innerRetina,os);

% irPlot(innerRetina, 'linear');
% irPlot(innerRetina, 'psth');

%% Look at covariance matrix
% load('/Users/james/Documents/MATLAB/isetbio misc/optimal linear decoder/ws_pixiumWhiteNoise_May4.mat')
psthstruct = mosaicGet(innerRetina.mosaic{1},'responsePsth');

psth1 = psthstruct.psth;
spikesout = psthstruct.spikes;

szCells = size(psth1);
cellCtr = 0;
szEnd = length(psth1{1,1});

whiteNoiseSmall = uint8(squeeze(whiteNoise.sceneRGB(:,:,:,1)));
% responseSpikes = mosaicGet(innerRetina.mosaic{1},'responseSpikes');   

spikesoutsm = uint8(spikesout);
% filename1 = ['/Users/james/Documents/MATLAB/isetbio misc/optimal linear decoder/May25_on2/WNstim_response_OnParasol_block_may25_' num2str(blockNum) '.mat'];
% filename1 = ['C:\Users\James\Documents\GitHub\offMidget1\WNstim_response_OffMidget_block_' num2str(blockNum) '.mat'];
filename1 = ['C:\Users\James\Documents\GitHub\onParasol1_all\WNstim_response_OnParasol_block_' num2str(blockNum) '.mat'];
% save(filename1, 'whiteNoiseSmall','spikesoutsm');
save(filename1, 'spikesoutsm');
toc
close
end