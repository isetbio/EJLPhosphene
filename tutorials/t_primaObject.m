% t_primaObject
% 
% The Prima object allows for the simulation of a subretinal prosthesis.
% The activation of the photovoltaic electrode array in response to a movie
% stimulus is generated. The bipolar current and RGC spikes are generated
% in response to the electrode stimulation. 

% With the RGC spikes, and a retina with a fixed set of parameters, we can
% generate the stimulus reconstruction to examine how it is affected by the
% prosthetic stimulation.

%% Load stimulus movie

stimSize = 100;

rsFactor = 1;

load([phospheneRootPath '/dat/stimuli/hallMovie.mat']);
szFrames = size(vidFrame,3);
hallMovieResize = zeros(rsFactor*stimSize,rsFactor*stimSize,szFrames);
for ii = 1:szFrames
    hallMovieResize(:,:,ii) = imresize(vidFrame(:,:,ii),[rsFactor*stimSize,rsFactor*stimSize]);
end

% movieIn = rand(100,100,10);
movieIn = hallMovieResize(:,:,1:100);

clear vidFrame hallMovieResize

%% Simulate bipolar and RGC response to prosthesis stimulation

primaParams.pixelWidth = 1*35e-6; % meters
primaParams.ecc = 1.8;       % mm
primaParams.fov = 1.7/1;     % deg

primaParams.pulseFreq = 100;           % Hz, electrode pulse frequency
primaParams.pulseDutyCycle = 1;        % Fraction of cycle pulse is on
primaParams.irradianceFraction = 1;    % Fraction of maximum irradiance 

primaRecon = primaArray(movieIn,primaParams)

primaRecon.compute(movieIn)

% Make a recon object here? Maybe primaRecon.reconstruct()?

%% Reconstruct - get spikes and decoding filter

spikeResp = mosaicSpikes(primaRecon.innerRetina);

rd = RdtClient('isetbio');
rd.crp('/resources/data/istim');
filterFile = 'filters_mosaic0_sv75_w1_sh2_may26primaSmall.mat';
data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
filterMat = data.filterMat; clear data;

%% Generate reconstructed movie
spikeAug(1,:) = ones(1,size(spikeResp,2));
spikeAug(2:9807,:) = spikeResp;
movRecon = filterMat'*spikeAug;
movReconPlay = reshape(movRecon,[100 100 size(spikeResp,2)]);
nFramesPlay = 40;
figure; ieMovie(movReconPlay(:,:,1:nFramesPlay));


%% Spatial zeroing of filter
% lambda = .01;
% filterMat2 = zeroFilter(filterMat,lambda);
% movRecon2 = filterMat2'*spikeAug;
% movReconPlay2 = reshape(movRecon2,[100 100 size(spikeResp,2)]);
% figure; ieMovie(movReconPlay2(:,:,1:nFramesPlay));
