% t_pixiumTestNS
% 
% Reconstructs the hallway navigation movie stimulus with decoder trained
% on the prosthesis simulation.
%
% Stimulating the activations of retinal ganglion cells with an array of
% electrodes. The stimulus image is defined, the electrode array is
% generated, the electrode activations are computed, the RGC mosaics are
% generated, the RGC mosaic responses are computed and the stimulus is
% inferred from the RGC mosaic resposnes using a linear decoder.
%
% Outline of computation:
% 1. Load image/movie
% 2. Outer segment representation
% 3. Build electrode array
% 4. Compute electrode activations from image/movie
% 5. Build RGC array
% 6. Calculate RGC input - only one, change to input from multiple
% 7. Build RGC activation functions
% 8. Compute RGC activations/spikes
% 9. Invert representation to form image/movie
% 10. Tile over big FOV
%
% 12/2016 JRG (c) isetbio
%
%  3bd0152 % nov 19\

%% Initialize
% clear;
% ieInit;
tic

%% Parameters to alter
clear electrodeArray

dropoutFlag = 0;

% Electrode size
% Set the size of implant pixels
electrodeArray.width = 70e-6; % meters
% electrodeArray.width = 140e-6; % meters

% % % % % % KNOBS
pulseFreq = 30;           % Hz, electrode pulse frequency
pulseDutyCycle = .3;       % Fraction of cycle pulse is on
irradianceFraction = .5;   % Fraction of maximum irradiance 
fadeTimeConstant = 5;      % Sets how quickly phosphene fades as f(energy); bigger time const means slower fade

percentDead = 0;

% %%% Grating subunit stimulus
%
% % params.barWidth = bwL;
% % iStim = ieStimulusGratingSubunit(params);
%
% % iStim = iStimC;
% absorptions = iStim.absorptions;
% % movingBar = iStim;
%
% nSteps = size(iStim.sceneRGB,3);
% params.nSteps = nSteps;
%
% size(iStim.sceneRGB)

% Retinal patch eccentricity
patchEccentricity = 1.5; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.6*2;

maxFreq = 30;

% Stimulus length
nSteps = 140;

params.nSteps = nSteps;
params.row = 100;
params.col = 100;
params.fov = fov;

%% Loop for different resizing factors
% The hall movies captures a large field of view. In order to reconstruct
% this large-FOV stimulus with a small RGC mosaic, we must tile the mosaic
% over the stimulus a number of times. The mosaic is necessarily small
% because the training algorithm takes a long time.

stimSize = 100;

for rsFactor = 1%[1 2 3 5 6]
    
    %% Resize the hallway movie stimulus for tiling
    rsFactor
    tic
    % Set full field movie stimulus

    numRows = 100; numCols = 100; gratingFreq = 4;
    gratingsMovie = buildGratings(numRows,numCols,gratingFreq);
    testmovieshort = (255*ieScale(gratingsMovie)); 
    
    % Stimulus parameters
    paramsStim.nsteps = 1;%nFrames;%size(testmovieshort,3);
    paramsStim.timeInterval = 1/125;%0.001; % sec
    paramsStim.expTime = 1/125;%0.001; % sec
    nFrames = nSteps;
    
    paramsStim.nsteps = nFrames;
    paramsStim.fov = 8;
    paramsStim.radius = 36/3*1e-6;
    paramsStim.theta = 330;
    paramsStim.side = 'left';
    paramsStim.fov = fov;
    blockctr = 0;
    
    tic
    
    %% Loop over each tiled mosaic
    for iblock = 1:rsFactor
        for jblock = 1:rsFactor
            
            % Get iStim structure with os object with desired properties
            blockctr = blockctr+1;
            paramsStim.nsteps = 1;
            iStim = ieStimulusMovie(testmovieshort((iblock-1)*stimSize+[1:stimSize],(jblock-1)*stimSize+[1:stimSize],1:10),paramsStim);
                        
            iStim.absorptions = iStim.sensor;
            clear movingBar
            movingBar = iStim;
            
            %% Outer segment calculation
            % There is no simulated outer segment, this identity outer segment acts as
            % a pass-through holder of the stimulus intensity information.
            
            % Input = RGB
            os = osCreate('displayrgb');
            
            % Get retinal patch properties for electrode array
            sceneSize = sceneGet(movingBar.scene,'size');
            retinalPatchWidth = sensorGet(movingBar.absorptions,'width','m');
            % retinalPatchWidth = sceneGet(movingBar.scene,'width');
            
            % retinalPatchHeight = sensorGet(movingBar.absorptions,'height','m');
            retinalPatchHeight = (sceneSize(1)/sceneSize(2))*retinalPatchWidth;
            
            % % % coneSpacing = scene.wAngular*300
            % coneSpacing = sensorGet(sensor,'dimension','um');
            os = osSet(os, 'patchSize', retinalPatchWidth);
            
            timeStep = sensorGet(movingBar.absorptions,'time interval','sec');
            os = osSet(os, 'timeStep', timeStep);
            
            movingBar.sceneRGB = testmovieshort((iblock-1)*stimSize+[1:stimSize],(jblock-1)*stimSize+[1:stimSize],:)
            % movingBar.sceneRGB = (contrastElectrode)*(movingBar.sceneRGB - 0.5)+0.5;
            os = osSet(os, 'rgbData', movingBar.sceneRGB);
            
            sceneRGB_Healthy = (1)*(movingBar.sceneRGB - 0.5)+0.5;
            osHealthy = os;
            osHealthy = osSet(osHealthy, 'rgbData', sceneRGB_Healthy);            
            
            retinalPatchSize = osGet(os,'size');
            
            % Electrode array properties
            metersPerPixel = retinalPatchWidth/retinalPatchSize(2);
            
            retinalPatchSize = osGet(os,'size');
            numberElectrodesX = floor(retinalPatchWidth/electrodeArray.width)+4;
            numberElectrodesY = floor(retinalPatchWidth/electrodeArray.width)+4;
            numberElectrodes = numberElectrodesX*numberElectrodesY;
            %% Build electrode array
            % Define the electrode array structure/object
            
            % Size stores the size of the array of electrodes
            electrodeArray.size = [numberElectrodesX numberElectrodesY];
            
            % Build the matrix of center coordinates for each electrode
            % electrodeArray.center(xPos,yPos,:) = [xCoord yCoord];
            x0 = -(numberElectrodesX-1)*electrodeArray.width/2;
            y0 = -(numberElectrodesY-1)*electrodeArray.width/2;
            for xPos = 1:numberElectrodesX
                for yPos = 1:numberElectrodesY
                    % electrodeArray.center(xPos,yPos,:) = [x0+(electrodeArray.width/2)*(xPos-1) + electrodeArray.width, y0+(electrodeArray.width/2)*(yPos-1) + electrodeArray.width];
                    electrodeArray.center(xPos,yPos,:) = [x0+(electrodeArray.width/1)*(xPos-1) + 0, y0+(electrodeArray.width/1)*(yPos-1) + 0 + (mod(xPos,2)-.5)*(electrodeArray.width/2)];
                end
            end
            
            th = (0:1/6:1)'*2*pi;
            xh = electrodeArray.width/2*cos(th);
            yh = electrodeArray.width/2*sin(th);
            
            % % Plot electrode array
            eaSize = size(electrodeArray.center);
            figure;
            hold on;
            for i = 1:eaSize(1)
                for j = 1:eaSize(2)
                    %         scatter(electrodeArray.center(i,j,1),electrodeArray.center(i,j,2));
                    plot(xh+electrodeArray.center(i,j,1),yh+electrodeArray.center(i,j,2),'r')
                end
            end
            axis equal
            xlabel('Distance (m)'); ylabel('Distance (m)');
            set(gca,'fontsize',14);
            
            % Build the current stimulation activation window
            % Gaussian activation from center of electrode
            activationWindow = floor(retinalPatchSize(2)/numberElectrodesX);
            electrodeArray.spatialWeight = fspecial('Gaussian', activationWindow, activationWindow/8);
            
            % Visualize Gaussian activation
            % figure; imagesc(electrodeArray.spatialWeight);
            % figure; surf(electrodeArray.spatialWeight);
            % xlabel(sprintf('Distance (\\mum)')); ylabel(sprintf('Distance (\\mum)'));
            % title('Gaussian Activation for a Single Electrode'); set(gca,'fontsize',16);
            
            %% Compute electrode activations from image
            
            % Get the full image/movie from the identity outersegment
            fullStimulus = osGet(os,'rgbData');
            
            % Find electrode activations by taking mean within window            
            for xPos = 1:numberElectrodesX
                for yPos = 1:numberElectrodesY
                    % Xcoords of window for stimulus
                    imageCoordX1 = (activationWindow)*(xPos-1)+1;
                    imageCoordX2 = (activationWindow)*(xPos);
                    
                    % Ycoords of window for stimulus
                    imageCoordY1 = (activationWindow)*(yPos-1)+1;
                    imageCoordY2 = (activationWindow)*(yPos);
                    
                    if imageCoordX2 > size(fullStimulus,2); imageCoordY2 = size(fullStimulus,2); end;
                    if imageCoordY2 > size(fullStimulus,1); imageCoordY2 = size(fullStimulus,1); end;
                    % Pull out piece of stimulus and take mean
                    electrodeStimulus = squeeze(fullStimulus(imageCoordY1:imageCoordY2,imageCoordX1:imageCoordX2,:,:));
                    electrodeArray.activation(xPos,yPos,:) = mean(RGB2XWFormat(electrodeStimulus));
                    
                    % sizeES = size(electrodeStimulus);
                    % electrodeArray.activation(xPos,yPos,frame) = min([ mean(electrodeStimulus(:,1:floor(sizeES(2)/2))) mean(electrodeStimulus(:,ceil(sizeES(2)/2):sizeES(2)))]);                    
                end
            end
                        
            %% Add X Hz spiking of stimulus
            
            % Right now the electrode sampling is at 0.01 s = 100 Hz
            % Downsample to get 5 Hz
            szAct = size(electrodeArray.activation);
            electrodeArray.activationDS = zeros(szAct);
            for iSample = 1:szAct(3)
                if mod(iSample,maxFreq/pulseFreq)==0
                    electrodeArray.activationDS(:,:,iSample) = electrodeArray.activation(:,:,iSample);
                    electrodeArray.activationDSoff(:,:,iSample) = 1-electrodeArray.activation(:,:,iSample);
                end
            end
            
            eaRS = reshape(electrodeArray.activation,[szAct(1)*szAct(2),szAct(3)]);
            eaDSRS = reshape(electrodeArray.activationDS,[szAct(1)*szAct(2),szAct(3)]);
            figure;
            % plot(eaRS');
            hold on;
            plot(eaDSRS');
            
            %% Build RGC array for healthy retina            
            clear paramsIR innerRetina
            if isunix || ismac
                %                 load([phospheneRootPath '/dat/mosaicAll_pix_ns.mat'])
                %                 load('mosaicAll_772530.mat');
                
                load('mosaicAll_35336498.mat');
%                 [a,b]=which('rdtclient');
%                 
%                 if isempty(a)
%                     load('mosaicAll_8372855.mat');
%                 else
%                     
%                     rdt = RdtClient('isetbio');
%                     rdt.crp('/resources/data/rgc');
%                     
%                     data = rdt.readArtifact('mosaicAll_8372855', 'type', 'mat');
%                     innerRetina = data.innerRetina;
%                 end
            else
%                 load([phospheneRootPath '\dat\mosaicAll_pix_ns.mat'])
            end
            %% Calculate RGC input
            % Weight electrode activation by Gaussian as a function of distance between
            % centers
            offFlag = [0 1 1 0]; nTileRows = 1; nTileCols = 1; mosaicOffset = zeros(4,1);
            innerRetinaInput = irActivationFromElectrode(innerRetina, electrodeArray, retinalPatchWidth, metersPerPixel, nTileRows, nTileCols, mosaicOffset, params, offFlag);
            
            %% Build RGC activation functions            
            innerRetinaThreshold = irGetThreshold(innerRetina);
            
            %% Compute RGC activations
            innerRetina =  irGetLinearRespElectrode(innerRetina, 1*innerRetinaInput, innerRetinaThreshold, params);
            
            %% Compute RGC spiking
            numberTrials = 1;
            for tr = 1:numberTrials
                innerRetina = irComputeSpikes(innerRetina,'coupling',false);
            end
            clear innerRetina1
            innerRetina1 = innerRetina;
            toc
            %% Do optimal reconstruction
            
            pOpt.innerRetina = innerRetina;
            pOpt.percentDead = percentDead;
            pOpt.numbins = 1;
%             pOpt.filterFile = 'pix1_long_filter_nsBig_100hz_4st';
%             pOpt.filterFile = 'pixiumBig/pix1_nsBig_100Hz_4st__mosaicAll_772530';
%             pOpt.filterFile = 'pixiumBig/pix1_nsBig_100Hz_1st_sv125__mosaicAll_772530';
%             pOpt.filterFile = 'pix1_nsBig_100Hz_1st_sv05__mosaicAll_8372855';
            pOpt.filterFile = 'pixium15_100/filters_pix1_nsBig15_1st_sv025__mosaicAll_35336498.mat';
       
            pOpt.percentDead = 0.8;
            [movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalReconSingle(pOpt);
            
            %% Save for tiling       
            
%             if ismac || isunix
%                 save([phospheneRootPath '/dat/ns_dec5_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetina','movrecons_on_offHealthy');
%             else
%                 save([phospheneRootPath '\dat\ns_dec5_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetina','movrecons_on_offHealthy');
%             end
            toc
        end
    end
    
    %% Tile reconstructed movies 
    szLen = size(testmovieshort,3);
    movieRecon = zeros(rsFactor*stimSize,rsFactor*stimSize,szLen);
    blockctr = 0; percentDead = 0; numbins = pOpt.numbins;
    tic
%     for iblock = 1:rsFactor
%         for jblock = 1:rsFactor
%             blockctr = blockctr+1;
%             if ismac || isunix
%                 load([phospheneRootPath '/dat/ns_dec5_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetina','movrecons_on_offHealthy');
%             else
% %                 load([phospheneRootPath '\dat\ns_dec5_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetina','movrecons_on_offHealthy');
%             end
% %             clear movrecons_on_offHealthy            
% %             pOpt.innerRetina = innerRetina;
% %             [movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalReconSingle(pOpt);
%             szLen = size(movrecons_on_offHealthy,3);
%             movieTmp =  movrecons_on_offHealthy;
%             movieRecon((iblock-1)*stimSize+[1:stimSize],(jblock-1)*stimSize+[1:stimSize],1:szLen) = movieTmp;% 255*ieScale(movieTmp - mean(movieTmp(:)));
%             clear movrecons_on_offHealthy
%         end
%     end
     
%     figure; ieMovie(movieRecon);
    
    %% Save as output with original movie
%     fig=figure;
%     opengl('software');
%     % set(fig,'position',[    624         437        1018         541]);
% %     set(fig,'position',[624   704   537   274]);
%     if ismac || isunix
%         aviobj = avifile([phospheneRootPath '/dat/prosthesis_dec20_recon_' num2str(rsFactor) '_ns.avi'])
%     else
%         aviobj = avifile([phospheneRootPath '\dat\prosthesis_dec20_recon_' num2str(rsFactor) '_ns.avi'])
%     end
%     aviobj.Fps = 30;
%     movieRecon = movrecons_on_offHealthy;

    integratedEnergy = [0:size(testmovieshort,3)-1]*irradianceFraction*pulseDutyCycle*(pulseFreq/maxFreq);
    fadeVal = ones(size(exp(-integratedEnergy/fadeTimeConstant)));

    shiftval = 1;
    szLen = size(testmovieshort,3)-shiftval;
    switch dropoutFlag
        case 0
            movieReconXW = RGB2XWFormat(movrecons_on_offHealthy(:,:,1:szLen+shiftval)).*(ones(size(movrecons_on_offHealthy,1)*size(movrecons_on_offHealthy,2),1)*fadeVal);
        case 1
            movieReconXW = RGB2XWFormat(movrecons_on_offHealthy_dropout(:,:,1:szLen+shiftval)).*(ones(size(movrecons_on_offHealthy,1)*size(movrecons_on_offHealthy,2),1)*fadeVal);
    end
    movieRecon = XW2RGBFormat(movieReconXW,size(movrecons_on_offHealthy,1),size(movrecons_on_offHealthy,2));
    movieComb = 255*irradianceFraction*pulseDutyCycle*ieScale(movieRecon(:,:,1:szLen-shiftval+1));
    movieComb(:,rsFactor*stimSize+[1:rsFactor*stimSize],:) = 255*irradianceFraction*pulseDutyCycle*ieScale(testmovieshort(:,:,shiftval+1:szLen+1));
%     movieComb(:,rsFactor*stimSize+[1:rsFactor*stimSize],:) = 255*ieScale(testmovieshort(:,:,shiftval+1:szLen+1));
    maxc = max(movieComb(:)); minc = min(movieComb(:));
%     for k=1:size(movieRecon,3)-60
%         % imagesc(movieRecon(:,:,k)); colormap gray;
%         
%         imagesc(movieComb(:,:,k)); colormap gray; axis image
%         
%         caxis([minc maxc]);
%         drawnow; pause(.001);
%         F = getframe(fig);
%         aviobj = addframe(aviobj,F);
%     end
%     close(fig)
%     aviobj = close(aviobj);
%     movieComb(1,1,:) = 255;
    %%
    figure;
    p.vname = [phospheneRootPath '/dat/pixiumBig/prosthesis_fullField_' num2str(rsFactor) '_ns_duty_' num2str(pulseDutyCycle*10) '_freq_' num2str(pulseFreq) '.avi']
    p.save = false;
    p.FrameRate = 30;
    ieMovie(movieComb, p);
    
    %%
    % clear movieComb movieRecon testmovieshort vidFrame
    toc
    
    %%
    
%     m1 = ieScale(movieRecon(:,:,1:567-shiftval+1));
%     m2 = ieScale(testmovieshort(:,:,shiftval+1:567+1));
%     err1 = sqrt((m1(:) - m2(:)).^2);
%     figure; hist(err1(:),40);
%     mean(err1(:)) %.170
%     std(err1(:)) % .141
    
end