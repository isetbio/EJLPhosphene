function movieOut = loadHallStimulus(stimFrames)


%% Load stimulus movie

stimSize = 100;

rsFactor = 1;

load(fullfile(phospheneRootPath, 'dat','stimuli','hallMovie.mat'));
szFrames = size(vidFrame,3);
hallMovieResize = zeros(rsFactor*stimSize,rsFactor*stimSize,szFrames);
for ii = 1:szFrames
    hallMovieResize(:,:,ii) = imresize(vidFrame(:,:,ii),[rsFactor*stimSize,rsFactor*stimSize]);
end

% movieIn = rand(100,100,10);
movieOut = hallMovieResize(:,:,1:stimFrames);

clear vidFrame hallMovieResize