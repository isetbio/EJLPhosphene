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

stimFrames = 50;
movieIn = loadHallMovie(stimFrames);

%% Simulate bipolar and RGC response to prosthesis stimulation

primaParams.pixelWidth = 1*35e-6; % meters
primaParams.ecc = 1.8;       % deg
primaParams.fov = 1.7/1;     % deg

primaParams.pulseFreq = 100;           % Hz, electrode pulse frequency
primaParams.pulseDutyCycle = 1;        % Fraction of cycle pulse is on
primaParams.irradianceFraction = 1;    % Fraction of maximum irradiance 

primaRecon = primaArray(movieIn,primaParams);

primaRecon.compute(movieIn)

%% Reconstruct - get spikes and decoding filter
 
spikeResp = mosaicSpikes(primaRecon.innerRetina);
% save('spikeResp_hallway.mat','spikeResp');

% Remote data toolbox - download decoding filter
rd = RdtClient('isetbio');
rd.crp('/resources/data/istim');
% filterFile = 'filters_mosaic0_sv75_w1_sh2_may26primaSmall';
% filterFile = 'filters_mosaic0_sv20_w1_sh2_dr0';
filterFile = 'filters_mosaic0_sv10_w1_sh2_dr0';
data  = rd.readArtifact(filterFile, 'type', 'mat');
filterMat = data.filterMat; clear data;

% % %% Generate reconstructed movie
spikeAug(1,:) = ones(1,size(spikeResp,2));
spikeAug(1+[1:size(spikeResp,1)],:) = spikeResp;
movRecon = filterMat'*spikeAug;
movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
nFramesPlay = stimFrames;
figure; ieMovie(movReconPlay(:,:,1:nFramesPlay));


% %% Spatial zeroing of filter
% lambda = .01;
% filterMat2 = zeroFilter(filterMat,lambda);
% movRecon2 = filterMat2'*spikeAug;
% movReconPlay2 = reshape(movRecon2,[100 100 size(spikeResp,2)]);
% figure; ieMovie(movReconPlay2(:,:,1:nFramesPlay));
