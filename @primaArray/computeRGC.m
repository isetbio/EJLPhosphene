function primaArray = computeRGC(primaArray)
%COMPUTERGC - a function of primaArray for computing the rgc layer spikes.
% 
%   primaArray.computeRGC();
% 
% At this point in the computation, we already have bipolar current due to
% electrode activations, so we can run the normal isetbio computation on
% RGC spikes from the bipolar mosaic.

%% Initialize the RGC mosaics

clear rgcL rgcParams

% Create retina ganglion cell layer object
rgcL = rgcLayer(primaArray.bpMosaic);

% There are various parameters you could set.  We will write a script
% illustrating these later.  We need a description.
rgcParams.centerNoise = 0;
rgcParams.ellipseParams = [1 1 0];  % Principle, minor and theta
% mosaicParams.axisVariance = .1;

% diameters = [3 3 1 1 5];  % In microns.

% 27*31+31*35+54*62+63*72
onPdiameter = 9.4;
diameters = [onPdiameter onPdiameter*.9 onPdiameter*.5 onPdiameter*.45];  % In microns.
    
cellType = {'on parasol','off parasol','on midget','off midget'};%,'onsbc'};
for ii = 1:length(cellType)
    rgcParams.rfDiameter = diameters(ii);
    rgcL.mosaic{ii} = rgcGLM(rgcL, primaArray.bpMosaic.mosaic{ii},cellType{ii},rgcParams);
end

nTrials = 1; rgcL.set('numberTrials',nTrials);

%% Compute the inner retina response and visualize

% Every mosaic has its input and properties assigned so we should be able
% to just run through all of them.
rgcL.compute('bipolarScale',250,'bipolarContrast',1);

primaArray.innerRetina = rgcL;
