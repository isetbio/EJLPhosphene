function initialize(primaArray, movieInput, varargin)
%INITIALIZE - a function of primaArray for initializing the object.


%%

p = inputParser;
addRequired(p, 'movieInput');
addParameter(p,  'bpFile','',@ischar);
addParameter(p,  'mosaicFile','',@ischar);
addParameter(p,  'fov',0,@isscalar);
addParameter(p,  'pixelWidth',0,@isscalar);
addParameter(p,  'height',0,@isscalar);
addParameter(p,  'width',0,@isscalar);
addParameter(p,  'ecc',0,@isscalar);
addParameter(p,  'pulseFreq',0,@isscalar);
p.KeepUnmatched = true;

p.parse(movieInput,varargin{:});
movieInput = p.Results.movieInput;
bpFile     = p.Results.bpFile;
mosaicFile = p.Results.mosaicFile;
fov        = p.Results.fov;
pixelWidth = p.Results.pixelWidth;
height      = p.Results.height;
width       = p.Results.width;
pulseFreq   = p.Results.pulseFreq;
ecc         = p.Results.ecc;

%% Cone Mosaic
% Generate the cone mosaic for retinal properties
% The absorptions and photocurrent will not be used

coneParams.fov = fov;
primaArray.fov = fov;

iStimNS = ieStimulusMovieCMosaic(rand(100,100,1),coneParams);
cMosaicNS = iStimNS.cMosaic;

%% Prima


primaArray.width = cMosaicNS.width;
primaArray.height = cMosaicNS.height;

primaArray.pixelWidth = pixelWidth;

primaArray.ecc = ecc;

primaArray.pulseFreq = pulseFreq;

retinalPatchSize = [cMosaicNS.height cMosaicNS.width];
retinalPatchHeight = cMosaicNS.height;
retinalPatchWidth = cMosaicNS.width;
% Electrode array properties
metersPerPixel = retinalPatchWidth/retinalPatchSize(2);

numberElectrodesX = floor(retinalPatchWidth/primaArray.pixelWidth)+4;
numberElectrodesY = floor(retinalPatchWidth/primaArray.pixelWidth)+4;
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

th = (0:1/6:1)'*2*pi;
xh = primaArray.pixelWidth/2*cos(th);
yh = primaArray.pixelWidth/2*sin(th);

% % Plot electrode array
eaSize = size(primaArray.center);
figure;
hold on;
for i = 1:eaSize(1)
    for j = 1:eaSize(2)
        %         scatter(primaArray.center(i,j,1),primaArray.center(i,j,2));
        plot(xh+primaArray.center(i,j,1),yh+primaArray.center(i,j,2),'r')
    end
end
axis equal
xlabel('Distance (m)'); ylabel('Distance (m)');
set(gca,'fontsize',14);

% Build the current stimulation activation window
% Gaussian activation from center of electrode
% activationWindow = floor(1e6*retinalPatchSize(2)/numberElectrodesX);
activationWindow = ceil(size(movieInput,1)/numberElectrodesX);
primaArray.spatialWeight = fspecial('Gaussian', round(1.5*activationWindow), 1.5*activationWindow/3);

% Visualize Gaussian activation
% figure; imagesc(primaArray.spatialWeight);
% figure; surf(primaArray.spatialWeight);
% xlabel(sprintf('Distance (\\mum)')); ylabel(sprintf('Distance (\\mum)'));
% title('Gaussian Activation for a Single Electrode'); set(gca,'fontsize',16);

%% Compute electrode activations from image
% 
% % Get the full image/movie from the identity outersegment
% fullStimulus = movieInput;
% 
% % Find electrode activations by taking mean within window
% for xPos = 1:numberElectrodesX
%     for yPos = 1:numberElectrodesY
%         % Xcoords of window for stimulus
%         imageCoordX1 = (activationWindow)*(xPos-1)+1;
%         imageCoordX2 = (activationWindow)*(xPos);
%         
%         % Ycoords of window for stimulus
%         imageCoordY1 = (activationWindow)*(yPos-1)+1;
%         imageCoordY2 = (activationWindow)*(yPos);
%         
%         if imageCoordX2 > size(fullStimulus,2); imageCoordY2 = size(fullStimulus,2); end;
%         if imageCoordY2 > size(fullStimulus,1); imageCoordY2 = size(fullStimulus,1); end;
%         % Pull out piece of stimulus and take mean
%         electrodeStimulus = squeeze(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2,:,:));
%         % primaArray.activation(xPos,yPos,:) = mean(RGB2XWFormat(electrodeStimulus));
%         
%         % Implement the local electrode min([e1,e2]) nonlinearity
%         sizeES = size(electrodeStimulus);
%         electrodeStimulusL = squeeze(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:floor(imageCoordX1+activationWindow/2),:,:));
%         electrodeStimulusR = squeeze(fullStimulus(imageCoordY1:imageCoordY2,floor(imageCoordX1+activationWindow/2)+1:imageCoordX2,:,:));
%         % primaArray.activation(xPos,yPos,frame) = min([ mean(electrodeStimulus(:,1:floor(sizeES(2)/2))) mean(electrodeStimulus(:,ceil(sizeES(2)/2):sizeES(2)))]);
%         primaArray.activation(xPos,yPos,:) = min([mean(RGB2XWFormat(electrodeStimulusL)); mean(RGB2XWFormat(electrodeStimulusL))]);
%     end
% end

%% Add X Hz spiking of stimulus
% primaArray.pulseFreq = pulseFreq;
% % Right now the electrode sampling is at 0.01 s = 100 Hz
% % Downsample to get 5 Hz
% szAct = size(primaArray.activation);
% primaArray.activationDS = zeros(szAct);
% for iSample = 1:szAct(3)
%     if mod(iSample,100/primaArray.pulseFreq)==0
%         primaArray.activationDS(:,:,iSample) = primaArray.activation(:,:,iSample);
%         primaArray.activationDSoff(:,:,iSample) = 1-primaArray.activation(:,:,iSample);
%     end
% end
% 
% eaRS = reshape(primaArray.activation,[szAct(1)*szAct(2),szAct(3)]);
% eaDSRS = reshape(primaArray.activationDS,[szAct(1)*szAct(2),szAct(3)]);

%% Bipolar
clear bpMosaic

cellType = {'ondiffuse','offdiffuse','onmidget','offmidget','onSBC'};
for cellTypeInd = 1:4
    clear bpParams
    bpParams.cellType = cellType{cellTypeInd};
    
    bpParams.ecc = primaArray.ecc;
    bpParams.rectifyType = 1;
    bpMosaic{cellTypeInd} = bipolar(cMosaicNS, bpParams);
    bpMosaic{cellTypeInd}.set('sRFcenter',1);
    bpMosaic{cellTypeInd}.set('sRFsurround',0);
%     bpMosaic{cellTypeInd}.compute(cMosaicNS);
end

primaArray.bpMosaic = bpMosaic;

%% RGC
clear params rgcParams
params.eyeRadius = primaArray.ecc;
params.eyeAngle = 90;
innerRetina=ir(bpMosaic,params);
cellType = {'on parasol','off parasol','on midget','off midget'};

rgcParams.centerNoise = 0;
rgcParams.model = 'LNP';
%     rgcParams.ellipseParams = [1 1 0];  % Principle, minor and theta

rgcParams.type = cellType{1};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{2};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{3};
innerRetina.mosaicCreate(rgcParams);
rgcParams.type = cellType{4};
innerRetina.mosaicCreate(rgcParams);

% innerRetina.compute(bpMosaic);

primaArray.innerRetina = innerRetina;
