function compute(primaArray, movieInput, varargin)
%COMPUTE - a function of primaArray for computing the object.
%       primaArray.compute(inputMovie);
% 
% This function computes the entire Prima pipeline:
%   - Movie stimulus input
%   - Prima electrode array activation
%   - Prima temporal sampling
%   - Bipolar stimulation by prima electrode array
%   - RGC stimulation by bipolar mosaic
% 
% See t_primaObject.m for a demonstration.
% 
% 5/2017 JRG (c) isetbio team
%% Movie input

p = inputParser;
addRequired(p, 'movieInput');
p.KeepUnmatched = true;

p.parse(movieInput,varargin{:});
movieInput = p.Results.movieInput;

%% Cone Mosaic
% Generate the dummy cone mosaic to get the properties of the retinal patch
% (size, etc.). The absorptions and photocurrent are not computed or used.

coneParams.fov = primaArray.fov;
iStimNS = ieStimulusMovieCMosaic(rand(100,100,1),coneParams);
cMosaicNS = iStimNS.cMosaic;

%% Compute electrode activations

primaArray.computeElectrode(movieInput);

%% Compute bipolar activations

primaArray.computeBipolar(cMosaicNS);

%% Compute RGC activations

primaArray.computeRGC();