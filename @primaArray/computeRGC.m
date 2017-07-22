function primaArray = computeRGC(primaArray)
%COMPUTERGC - a function of primaArray for computing the rgc layer spikes.
% 
%   primaArray.computeRGC();
% 
% At this point in the computation, we already have bipolar current due to
% electrode activations, so we can run the normal isetbio computation on
% RGC spikes from the bipolar mosaic.

%% Initialize the RGC mosaics
clear params rgcParams
params.eyeRadius = primaArray.ecc;
params.eyeAngle = 90;
innerRetina=ir(primaArray.bpMosaic,params);
cellType = {'on parasol','off parasol','on midget','off midget'};

rgcParams.centerNoise = 0;
rgcParams.model = 'LNP';
% rgcParams.ellipseParams = [1 1 0];  % Principle, minor and theta
% rng(20001);
rgcParams.axisVariance = 0;%.085;
rgcParams.centerNoise = 0;%.05;
rgcParams.type = cellType{1};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{2};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{3};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{4};
innerRetina.mosaicCreate(rgcParams);

scaleFactor = 1;
innerRetina = scaleRF(innerRetina, scaleFactor);

% Compute spikes

innerRetina.compute(primaArray.bpMosaic);

primaArray.innerRetina = innerRetina;
