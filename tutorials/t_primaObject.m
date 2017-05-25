% t_primaObject

% Generate stimulus
% Build primaArray
%   pseudo cMosaic
%   electrode array
%   electrode activation
%   electrode threshold
%   bipolar activation
%   rgc spikes
% 
% return bp, ir

% do recon


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
movieIn = hallMovieResize(:,:,1:250);
%%

primaParams.pixelWidth = 1*35e-6; % meters
primaParams.ecc = 1.8;       % mm
primaParams.fov = 1.7/1;     % deg

primaParams.pulseFreq = 100;           % Hz, electrode pulse frequency
primaParams.pulseDutyCycle = 1;        % Fraction of cycle pulse is on
primaParams.irradianceFraction = 1;    % Fraction of maximum irradiance 

primaRecon = primaArray(movieIn,primaParams)

primaRecon.compute(movieIn)

% primaRecon.computeBipolar()
% 
% primaRecon.computeRGC()


