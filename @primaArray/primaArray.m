classdef primaArray < hiddenHandle
%PRIMAARRAY - Create a primaArray object
% The primaArray class is used to simulate the Pixium Prima subretinal
% prosthesis.
% 
% primaPros = primaArray(isetbioObject, 'PARAM1', val1, 'PARAM2', val2,...)
% creates the primaArray object. Optional parameter name/value pairs are
% listed below.
% 
% The primaArray object computes the response of bipolar cells and RGCs
% stimulated by the Prima prosthesis.
%
%       pPrima.ecc = 1;
%       primaPros = primaArray(pRecon);
% 
% Input: 
% 
% Output: the inner retina
% 
%     cellLocation;                    % location of bipolar RF center
% 
%  ISETBIO wiki: <a href="matlab:
%  web('https://github.com/isetbio/isetbio/wiki/prima','-browser')">reconstruction</a>.
%   
% 5/2017 JRG (c) isetbio team

%% Define object
% Public, read-only properties.
properties (SetAccess = private, GetAccess = public)
end

% Protected properties.
properties (SetAccess = protected, GetAccess = public)
end

% Private properties. Only methods of the parent class can set these
properties(Access = private)
end

properties (Access = public)    
    %MOSAICFILE - name of file used to generate inner retina mosaic
    mosaicFile;
    
    %SIZE - in number of electrodes
    size;
    
    %PIXELWIDTH - width of pixel
    pixelWidth;
    
    %WIDTH - width of array in meters
    width;
    
    %HEIGHT - height of array in meters
    height;
    
    %CENTER - array of center locations
    center;
    
    %FOV - field of view in degrees of visual angle
    fov;
    
    %ECC - eccentricity in degrees of visual angle
    ecc;
    
    %SPATIALWEIGHT - activation weight
    spatialWeight;
    
    %PULSEFREQ - frequency of electrode activation
    pulseFreq;
       
    %PULSEDUTYCYCLE - fraction of cycle pulse is on
    pulseDutyCycle;        
    
    %IRRADIANCEFRACTION - fraction of maximum irradiance
    irradianceFraction;    
    
    %CURRENTDECAY - stdev of gaussian that weighs current spread from electrode
    currentDecay;
    
    %ACTIVATION - electrode activations in arbitray units
    activation;
    activationDS;
    activationDSoff;
    
    %BPMOSAIC - holds bipolar objects
    bpMosaic;
    
    %INNERRETINA - holds inner retina object
    innerRetina;
    
    electrodeCoords;
    electrodeCoordsFull;
end

% Public methods
methods
    
    % Constructor
    function obj = primaArray(movieInput,varargin)     
        % Initialize the primaArray class
        %   primaPros = primaArray();
        
        p = inputParser;
        addRequired(p, 'movieInput');
        addParameter(p,  'bpFile','',@ischar);
        addParameter(p,  'mosaicFile','',@ischar);
        
        p.KeepUnmatched = true;
        
        p.parse(movieInput,varargin{:});  
        
        bpFile = p.Results.bpFile;
        mosaicFile = p.Results.mosaicFile;
                      
        obj.initialize(movieInput,p.Unmatched)
    end
    
    obj = initialize(obj, movieInput, primaParams);
    
    % Declare the method for calculating electrode activation
    computeElectrode(primaArray, movieInput);
    
    computeBipolar(primaArray, cMosaic);    
    
    computeRGC(primaArray);    
           
    visualizeStimulusAndElectrodeActivation(primaArray, filename, fullStimulus, linearActivation, activation, activationDS, activationDSoff);
    visualizePhotocurrentAndBpMosaicResponses(primaArray, filename, weights, photocurrentResponse, bpResponseCenter, bpResponseCenterFull);

%     function window(obj)
%         obj.figureHandle = primaWindow(obj);
%         % Tip: Retrieve guidata using
%         %    gui = guidata(obj.figureHandle);
%         %
%     end
    
end

properties (Constant)
end

 % Methods may be called by the subclasses, but are otherwise private
methods (Access = protected)
end

% Methods that are totally private (subclasses cannot call these)
methods (Access = private)
end

end