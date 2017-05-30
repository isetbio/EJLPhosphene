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

% Electrode size
% Set the size of implant pixels
electrodeArray.width = 25e-6; % meters
% electrodeArray.width = 140e-6; % meters

% Retinal patch eccentricity
patchEccentricity = 4; % mm

% Field of view/stimulus size
% Set horizontal field of view
fov = 1.6*2;

% % % % % % KNOBS
pulseFreq = 100;           % Hz, electrode pulse frequency
pulseDutyCycle = .2;       % Fraction of cycle pulse is on
irradianceFraction = .5;  % Fraction of maximum irradiance 

% Stimulus length
nSteps = 520;

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
%     load('C:\Users\James\Documents\MATLAB\github\EJLPhosphene\dat\stimuli\hallMovie.mat')
    load([phospheneRootPath '/dat/stimuli/hallMovie.mat'])
    szFrames = size(vidFrame,3);
    hallMovieResize = zeros(rsFactor*stimSize,rsFactor*stimSize,szFrames);
    for ii = 1:szFrames
        hallMovieResize(:,:,ii) = imresize(vidFrame(:,:,ii),[rsFactor*stimSize,rsFactor*stimSize]);
    end
    
    % Set hallway movie stimulus
    testmovieshort = (255*ieScale(hallMovieResize)); clear hallMovieResize;
    
    % Stimulus parameters
    paramsStim.nsteps = 1;%nFrames;%size(testmovieshort,3);
    paramsStim.timeInterval = 1/125;%0.001; % sec
    paramsStim.expTime = 1/125;%0.001; % sec
    nFrames = nSteps;
    
    paramsStim.nsteps = nFrames;
    paramsStim.fov = 2*1.6;
    paramsStim.radius = 36/3*1e-6;
    paramsStim.theta = 330;
    paramsStim.side = 'left';
    paramsStim.fov = fov;
    blockctr = 0;
    
    tic
    
    %% Loop over each tiled mosaic
    iblock = 1;
        jblock = 1;
            
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
            
            % movingBar.sceneRGB = testmovieshort((iblock-1)*stimSize+[1:stimSize],(jblock-1)*stimSize+[1:stimSize],:)
            % movingBar.sceneRGB = (contrastElectrode)*(movingBar.sceneRGB - 0.5)+0.5;
%             
%             movingBar.sceneRGB = repmat(imgLandoltC,[1 1 7]);
%             os = osSet(os, 'rgbData', movingBar.sceneRGB);
%             
%             sceneRGB_Healthy = (1)*(movingBar.sceneRGB - 0.5)+0.5;
%             osHealthy = os;
%             osHealthy = osSet(osHealthy, 'rgbData', sceneRGB_Healthy);            
            
        
            %% Build RGC array for healthy retina            
            clear paramsIR innerRetina
            if isunix || ismac
%                 load([phospheneRootPath '/dat/mosaicAll_pix_ns.mat'])
%                 load('mosaicAll_772530.mat');
%                 load('mosaicAll_8372855.mat');
%                     load('mosaicAll_19261.mat');
%                     load('mosaicAll_51550348.mat');
                    load('mosaicAll_1246640.mat');
%                     load('mosaicAll_20116.mat');
            else
%                 load([phospheneRootPath '\dat\mosaicAll_pix_ns.mat'])
            end
%%
            tic
            orient0 = 'down';
%     for expRep = 1:500;
        cFrames = 7;
%         expRep
            movingBar.sceneRGB = repmat(imgLandoltC('orientation',orient0),[1 1 cFrames]);
%             os = osSet(os, 'rgbData', movingBar.sceneRGB);

        contarr = [.01 .025 .05 .005];
        
        for contind = 1:length(contarr)

            sceneRGB_Healthy = (contarr(contind))*(movingBar.sceneRGB - 0.5);%+0.5;
            osHealthy = os;
            osHealthy = osSet(osHealthy, 'rgbData', sceneRGB_Healthy); 
    for expRep = 1:500;
        cFrames = 7;
        expRep
            innerRetina = irCompute(innerRetina,osHealthy,'coupling',false);
           
%             toc
            %% Do optimal reconstruction
            
            pOpt.innerRetina = innerRetina;
            pOpt.percentDead = 0;
%             pOpt.numbins = 4;
%             pOpt.filterFile = 'pix1_long_filter_nsBig_100hz_4st';
%             pOpt.filterFile = 'pixiumBig/pix1_nsBig_100Hz_4st__mosaicAll_772530';
%             pOpt.filterFile = 'pixiumBig/pix1_nsBig_100Hz_1st_sv125__mosaicAll_772530';
%             pOpt.filterFile = 'pixium25/pix1_nsBig_100Hz_1st_sv05__mosaicAll_8372855';


%             mosaicFile = '_mosaicAll_19261';
%             mosaicFile = '_mosaicAll_51550348';
            mosaicFile = '_mosaicAll_1246640';
%             mosaicFile = '_mosaicAll_20116';
% 
            pOpt.numbins = 1;
%             filterFile = ['ns100/filters_nsDec22_1st_sv05_' mosaicFile];
%                 filterFile = ['ns100_r2/filters_ns100_Dec31_1st_sv05_' mosaicFile];
                
%             filterFile = ['ns100_r2_10/filters_ns100_jan1_1st_sv20_' mosaicFile];
%             pOpt.numbins = 2;
%             filterFile = ['ns200/filters_nsDec22_2st_sv375fig_' mosaicFile];

            pOpt.numbins = 1;
%             filterFile = ['ns200/filters_nsDec22_4st_sv075_' mosaicFile];
%             filterFile = ['ns200/filters_nsDec22_4st_sv005_' mosaicFile];

%             filterFile = ['ns100_r2_10/filters_ns1002_jan1_1st3_sv20_' mosaicFile];
            
%             filterFile = ['ns100_r2_10/filters_ns100_jan1_1st_sv05_' mosaicFile];
            
%             filterFile = ['ns100_r2_10/filters_ns100_jan1_1st_sh9_sv20_' mosaicFile];
             filterFile = ['ns100_r2_10/filters2_ns100_feb6_sh9_sv40_tr83' mosaicFile];
            
%             filterFile = ['ns100_r2_10_regmos/filters_ns100_regmos_1st_sh9_sv30_' mosaicFile];
%             filterFile = ['ns100/filters_nsDec22_4st_sv05_' mosaicFile];
% 
%             pOpt.numbins = 8;
% % %             filterFile = ['ns200/filters_nsDec22_8st_sv125_' mosaicFile];
%             filterFile = ['ns200/filters_nsDec22_8st_sv03_' mosaicFile];

            pOpt.filterFile = filterFile;
            [movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalReconSingle(pOpt);
%             figure; ieMovie(movrecons_on_offHealthy);
            movrecons_on_offHealthy = movrecons_on_offHealthy(:,:,1:cFrames);
            %% Save for tiling       
            
%             if ismac || isunix
%                 save([phospheneRootPath '/dat/ns_dec5_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetina','movrecons_on_offHealthy');
%             else
%                 save([phospheneRootPath '\dat\ns_dec5_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetina','movrecons_on_offHealthy');
%             end
            toc
            
            framesOut(:,:,expRep) = movrecons_on_offHealthy(:,:,1);
            
%             if expRep == 250
%                  save(['dat/framesOutLong' orient0 '_250_c075_418.mat'],'framesOut');
%             end
                
    end
    toc
    save(['dat/framesOutLong' orient0 '_c0' num2str(1000*contarr(contind)) '_420.mat'],'framesOut');
        end
        toc
    
    
    %% Tile reconstructed movies 
    szLen = size(movrecons_on_offHealthy,3);
%     movieRecon = zeros(rsFactor*stimSize,rsFactor*stimSize,szLen);
%     blockctr = 0; percentDead = 0; numbins = pOpt.numbins;
%     tic
%     for iblock = 1:rsFactor
%         for jblock = 1:rsFactor
%             blockctr = blockctr+1;
%             if ismac || isunix
%                 load([phospheneRootPath '/dat/pixiumBig/ns_dec20_rs_' num2str(rsFactor) '_block' num2str(blockctr) '.mat'], 'innerRetina','movrecons_on_offHealthy');
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
%      
% %     figure; ieMovie(movieRecon);
    
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
% % %     movieRecon = movrecons_on_offHealthy;
% % %     shiftval = 9;
% % %     clear movieComb
% % %     szLen = 596;
% % %     movieComb = 255*irradianceFraction*pulseDutyCycle*ieScale(movieRecon(:,:,1:szLen-shiftval+1));
% % %     movieComb(:,rsFactor*stimSize+[1:rsFactor*stimSize],:) = 255*irradianceFraction*pulseDutyCycle*ieScale(testmovieshort(:,:,shiftval+1:szLen+1));
% % %     maxc = max(movieComb(:)); minc = min(movieComb(:));
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
    
    %%
%     figure;
%    % p.vname = [phospheneRootPath '/dat/pixiumBig/prosthesis_dec20_recon_' num2str(rsFactor) '_ns_fr25sum.avi']
%     p.vname = [reconstructionRootPath '/dat/ns100_r2/ns_Dec31_recon_' num2str(rsFactor) '_8st_ev125.avi']
%   
%     p.save = false;
%     p.FrameRate = 25;
%     ieMovie(movieComb, p);
    
    %%
    
%     mc1 = ieScale(movieRecon(:,:,1:szLen-shiftval+1));
%     mc2 = ieScale(testmovieshort(:,:,shiftval+1:szLen+1));
%     errmov =mc1-mc2;
%     errtot = ((errmov.^2));
%     figure; ieMovie(errmov);
%     figure; subplot(131); imagesc(mc1(:,:,1)); subplot(132); imagesc(mc2(:,:,1)); subplot(133); imagesc(mc1(:,:,1)-mc2(:,:,1));
    
%     3: rms = .156, rss = .1981
%       1: .1596, .202
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