function initialize(primaArray, movieInput, varargin)
%INITIALIZE - a function of primaArray for initializing the object.
% 
%   primaRecon = primaArray(movieIn,primaParams);
% 
% This function generates the primaArray object. The general approach is to
% simulate the electrode array, the bipolar layer and the rgc layer. Since
% we are simulating a subretinal prosthesis, we assume the cone
% photoreceptors are not functional. However, we create a dummy cone mosaic
% that is only used to simulate a retinal patch and to build the bipolar
% layer. The response of the cone mosaic to the stimulus is not computed.
% 
% The general outline of the initiailization is:
%   - Get input movie
%   - Build dummy cone mosaic to get retinal properties at ecc and fov
%   - Bulid the prima electrode array with properties from dummy cone mosaic
%   - Build the bipolar layer
%   - Build the rgc layer
% 
% Required paramters are:
%   - movieInput: the movie stimulus, MxNxK matrix
% 
% Optional parameters that may be set include:
%   - bpFile: a saved bipolar layer (in case trained recon filters had irregular mosaic)
%   - mosaicFile: a saved rgc layer (in case trained recon filters had irregular mosaic)
%   - fov: field of view in degrees
%   - pixelWidth: the pitch between pixels in the electrod array (microns)
%   - height: height of array
%   - width: width of array
%   - ecc: eccentricity of retinal patch (degrees)
%   - pulseFreq: the pulse rate of the electrode array (Hz)
%   - pulseDutyCycle: the duty cycle of electrode pulse (%)
%   - irradianceFraction: the intensity of light projected onto electrodes (%)
%   - currentDecay: the spatial falloff of current stimulation on bipolars 
% 
% 5/2017 JRG (c) isetbio team

%% Input parser

p = inputParser;
addRequired(p, 'movieInput');
addParameter(p,  'bpFile','',@ischar);
addParameter(p,  'mosaicFile','',@ischar);
addParameter(p,  'fov',0,@isscalar);
addParameter(p,  'pixelWidth',70,@isscalar);
addParameter(p,  'height',0,@isscalar);
addParameter(p,  'width',0,@isscalar);
addParameter(p,  'ecc',0,@isscalar);
addParameter(p,  'pulseFreq',0,@isscalar);
addParameter(p,  'pulseDutyCycle',1,@isscalar);
addParameter(p,  'irradianceFraction',1,@isscalar);
addParameter(p,  'currentDecay',2,@isscalar);
p.KeepUnmatched = true;

p.parse(movieInput,varargin{:});
movieInput = p.Results.movieInput;
bpFile     = p.Results.bpFile;
mosaicFile = p.Results.mosaicFile;
fov        = p.Results.fov;
pixelWidth = p.Results.pixelWidth;
height      = p.Results.height;
width       = p.Results.width;
ecc         = p.Results.ecc;
pulseFreq   = p.Results.pulseFreq;
pulseDutyCycle     = p.Results.pulseDutyCycle;
irradianceFraction = p.Results.irradianceFraction;
currentDecay      = p.Results.currentDecay;


%% Cone Mosaic
% Generate the dummy cone mosaic to get the properties of the retinal patch
% (size, etc.). The absorptions and photocurrent are not computed or used.

coneParams.fov = fov;
primaArray.fov = fov;

% We generate the cone mosaic with only one frame of the stimulus to save
% time on the computation.
iStimNS = ieStimulusMovieCMosaic(rand(100,100,1),coneParams);
cMosaicNS = iStimNS.cMosaic;

%% Prima
% We initialize many of the primaArray object properties from the
% primaParams input structure.

% From input
primaArray.pixelWidth = pixelWidth;
primaArray.ecc = ecc;
primaArray.pulseFreq = pulseFreq;
primaArray.pulseDutyCycle = pulseDutyCycle;
primaArray.irradianceFraction = irradianceFraction;
primaArray.currentDecay = currentDecay;

% From the dummy cone mosaic
primaArray.width = cMosaicNS.width;
primaArray.height = cMosaicNS.height;
retinalPatchSize = [cMosaicNS.height cMosaicNS.width];
retinalPatchHeight = cMosaicNS.height;
retinalPatchWidth = cMosaicNS.width;

% Electrode array properties derived from dummy cone mosaic
metersPerPixel = retinalPatchWidth/retinalPatchSize(2);

% Get size of electrode array
numberElectrodesX = floor(retinalPatchWidth/primaArray.pixelWidth)+1;
numberElectrodesY = floor(retinalPatchWidth/primaArray.pixelWidth)+1;
numberElectrodes = numberElectrodesX*numberElectrodesY;


%% Build electrode array
% Define the electrode array structure/object

% Size stores the size of the array of electrodes
primaArray.size = [numberElectrodesX numberElectrodesY];

% Build the matrix of center coordinates for each electrode
% primaArray.center(xPos,yPos,:) = [xCoord yCoord];
x0 = -(numberElectrodesX-1)*primaArray.pixelWidth/2;
y0 = -(numberElectrodesY-1)*primaArray.pixelWidth/2;
for xPos = 1:numberElectrodesX
    for yPos = 1:numberElectrodesY
        % primaArray.center(xPos,yPos,:) = [x0+(primaArray.pixelWidth/2)*(xPos-1) + primaArray.pixelWidth, y0+(primaArray.pixelWidth/2)*(yPos-1) + primaArray.pixelWidth];
        primaArray.center(xPos,yPos,:) = [x0+(primaArray.pixelWidth/1)*(xPos-1) + 0, y0+(primaArray.pixelWidth/1)*(yPos-1) + 0 + (mod(xPos,2)-.5)*(primaArray.pixelWidth/2)];
    end
end

%% Build the current stimulation activation window
% Gaussian activation from center of electrode
% activationWindow = floor(1e6*retinalPatchSize(2)/numberElectrodesX);
activationWindow = ceil(size(movieInput,1)/numberElectrodesX);
primaArray.spatialWeight = fspecial('Gaussian', round(1.5*activationWindow), 1.5*activationWindow/3);

%% Initialize the bipolar layer

clear bpMosaic

% Build each subtype of bipolar mosaic
cellType = {'ondiffuse','offdiffuse','onmidget','offmidget','onSBC'};
for cellTypeInd = 1:4
    clear bpParams
    bpParams.cellType = cellType{cellTypeInd};
    
    bpParams.ecc = primaArray.ecc;
    bpParams.rectifyType = 1;
    bpMosaic{cellTypeInd} = bipolar(cMosaicNS, bpParams);
    bpMosaic{cellTypeInd}.set('sRFcenter',1);
    bpMosaic{cellTypeInd}.set('sRFsurround',0);
end
% Attach to the prima object
primaArray.bpMosaic = bpMosaic;

%% Initialize the RGC layer

% RGC parameters
clear params rgcParams
params.eyeRadius = primaArray.ecc;
params.eyeAngle = 90;
innerRetina=ir(bpMosaic,params);
cellType = {'on parasol','off parasol','on midget','off midget'};

% Choose how regular the mosaic is
rgcParams.centerNoise = 0;
rgcParams.model = 'LNP';
%     rgcParams.ellipseParams = [1 1 0];  % Principle, minor and theta

% Generate each mosaic of rgc cell type
rgcParams.type = cellType{1};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{2};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{3};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{4};
innerRetina.mosaicCreate(rgcParams);

% Attach to the primaArray object
primaArray.innerRetina = innerRetina;

%% Make a figure for the electrode array
% th = (0:1/6:1)'*2*pi;
% scaleFactor = 1e6;
% xh = scaleFactor*primaArray.pixelWidth/2*cos(th);
% yh = scaleFactor*primaArray.pixelWidth/2*sin(th);
% 
% % % Plot electrode array
% eaSize = size(primaArray.center);
% % figure;
% hold on;
% for i = 1:eaSize(1)
%     for j = 1:eaSize(2)
%         %         scatter(primaArray.center(i,j,1),primaArray.center(i,j,2));
%         plot(xh+scaleFactor*primaArray.center(i,j,1),yh+scaleFactor*primaArray.center(i,j,2),'r')
%     end
% end
% axis equal
% xlabel('Distance (m)'); ylabel('Distance (m)');
% set(gca,'fontsize',14);
