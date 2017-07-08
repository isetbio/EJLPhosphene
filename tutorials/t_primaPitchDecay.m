% t_primaObject
% 
% The Prima object allows for the simulation of a subretinal prosthesis.
% The activation of the photovoltaic electrode array in response to a movie
% stimulus is generated. The bipolar current and RGC spikes are generated
% in response to the electrode stimulation. 
% 
% With the RGC spikes, and a retina with a fixed set of parameters, we can
% generate the stimulus reconstruction to examine how it is affected by the
% prosthetic stimulation.
% 
% 
% TOOLBOX DEPENDENCIES - these must be downloaded and added to the matlab
%                               path with subfolders.
%       isetbio:            http://github.com/isetbio/isetbio [bipolar branch]
%       RGC-Reconstruction: https://github.com/Chichilnisky-Lab/RGC-Reconstruction
%       EJLPhosphene:       https://github.com/isetbio/EJLPhosphene
%       RemoteDataToolbox:  https://github.com/isetbio/RemoteDataToolbox
% 

%% Load stimulus movie

stimFrames = 550;
movieIn = loadHallStimulus(stimFrames);

%% Simulate bipolar and RGC response to prosthesis stimulation

% primaParams.pixelWidth = 1*35e-6; % meters
primaParams.ecc = 1.8;       % deg
primaParams.fov = 1.7/1;     % deg

primaParams.pulseFreq = 100;           % Hz, electrode pulse frequency
primaParams.pulseDutyCycle = 1;        % Fraction of cycle pulse is on
primaParams.irradianceFraction = 1;    % Fraction of maximum irradiance 

% 
pixelPitchArr = [70 70/2 70/4 70/8]*1e-6;% 70/16 70/32];
currentDecayArr = [2 4 8 16];

for pitchInd = [1:length(pixelPitchArr)]
    for currentDecayInd = 1:length(currentDecayArr)
[pitchInd currentDecayInd]
primaParams.pixelWidth = pixelPitchArr(pitchInd);
primaParams.currentDecay = currentDecayArr(currentDecayInd);

primaRecon = primaArray(movieIn,primaParams);

primaRecon.compute(movieIn)

%% Reconstruct - get spikes and decoding filter
 
spikeResp = mosaicSpikes(primaRecon.innerRetina);
save(['spikeResp_hallway_pitch' num2str(pitchInd) '_decay' num2str(currentDecayInd) '.mat'],'spikeResp');

    end
end


% Remote data toolbox - download decoding filter
% rd = RdtClient('isetbio');
% rd.crp('/resources/data/istim');
% % filterFile = 'filters_mosaic0_sv75_w1_sh2_may26primaSmall';
% filterFile = 'filters_mosaic0_sv20_w1_sh2_dr0';
% % filterFile = 'filters_mosaic0_sv10_w1_sh2_dr0';
% data  = rd.readArtifact(filterFile, 'type', 'mat');
% filterMat = data.filterMat; clear data;

filterFile = 'filters_mosaic0_sv20_w1_sh2_dr0';
load(filterFile);

% % %% Generate reconstructed movie
spikeAug(1,:) = ones(1,size(spikeResp,2));
spikeAug(1+[1:size(spikeResp,1)],:) = spikeResp;
movRecon = filterMat'*spikeAug;
movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
nFramesPlay = stimFrames;
% figure; ieMovie(movReconPlay(:,:,1:nFramesPlay));


%% Spatial zeroing of filter
% lambda = .025;
% filterMat2 = zeroFilter(filterMat,lambda);
% movRecon2 = filterMat2'*spikeAug;
% movReconPlay2 = reshape(movRecon2,[100 100 size(spikeResp,2)]);
% figure; ieMovie(movReconPlay2(:,:,1:nFramesPlay));

%% No Learning
%  27*31+31*35+54*62+63*72
onParasolInd = 1+27*31;
offParasolInd = onParasolInd+31*35;
onMidgetInd = offParasolInd+54*62;
load('filters_mosaic0_sv50_w1_sh4_tr80.mat')
spikeNoLearn = spikeAug;
spikeNoLearn(onParasolInd+1:offParasolInd,:) = 0;
spikeNoLearn(onMidgetInd+1:end,:) = 0;

% spikeNoLearn(1:onParasolInd,:) = 0;
% spikeNoLearn(offParasolInd+1:onMidgetInd,:) = 0;

movReconNoLearn = filterMat'*spikeNoLearn;
movReconPlayNoLearn = reshape(movReconNoLearn,[100 100 size(spikeResp,2)]);
movReconPlayNoLearn = permute(movReconPlayNoLearn,[2 1 3]);
nFramesPlay = stimFrames;
% figure; ieMovie(movReconPlayNoLearn(:,:,1:nFramesPlay));

%%
lambda = .01;
filterMat2 = zeroFilter(filterMat,lambda);
movReconNoLearn2 = filterMat2'*spikeNoLearn;
movReconPlayNoLearn2 = reshape(movReconNoLearn2,[100 100 size(spikeResp,2)]);

movReconPlayNoLearn2 = permute(movReconPlayNoLearn2,[2 1 3]);
% figure; ieMovie(movReconPlayNoLearn2(:,:,1:nFramesPlay));


