function [movrecons_on_off_full, movrecons_on_off_dropout, spikeRespOnOff] = irOptimalReconSingle(varargin)
% Converts an RGC mosaic spike train to a reconstructed movie.
% 
% testReconAll loads an RGC mosaic as well as a movie input and computes 
% the spike response. The spike response is passed to this function, where
% the spikes are pulled out of the object, rearranged and linearly decoded
% to produce a reconstructed movie.
% 
% 10/2016 JRG (c) isetbio team
% 
% See also testReconAll, runReconstruct, irOptimalRecon

%%


p = inputParser;
p.addParameter('innerRetina',[]);
p.addParameter('percentDead',0,@isnumeric);
p.addParameter('filterFile',[],@ischar);
p.addParameter('numbins',[],@isnumeric);
p.parse(varargin{:});
innerRetina = p.Results.innerRetina;
percentDead = p.Results.percentDead;
filterFile = p.Results.filterFile;
numbins = p.Results.numbins;
%% Get spikes from mosaic objects

y = cell(length(innerRetina.mosaic));

dt = innerRetina.mosaic{1}.dt;

% Loop over the mosaic objects

maxTrials = innerRetina.numberTrials;
for trial = 1:maxTrials
    for mosaicInd = 1:length(innerRetina.mosaic)
        
        cellCtr=0;
        nCells = size(innerRetina.mosaic{mosaicInd}.responseSpikes);
        
        y{mosaicInd} = zeros(nCells(1)*nCells(2),10000);
        
        for ycell = 1:nCells(1)
            for xcell = 1:nCells(2)
                cellCtr = cellCtr+1;
                yind =  innerRetina.mosaic{mosaicInd}.responseSpikes{xcell,ycell,trial,1};
                y{mosaicInd}(cellCtr,ceil(yind./dt))=y{mosaicInd}(cellCtr,ceil(yind./dt))+1;
            end
        end
    end
end


%% Load the reconstruction filters, locally or from RDT

% all four mosaics, nov 4
% load('/Users/james/Downloads/filters_mosaic_ns_all_overlap0_svd_1000_len_100.mat')
% % load('/Users/james/Downloads/filters_mosaic_wn_all_42reps_overlap0_svd_1000_len_100.mat')
% 
% load('C:\Users\James\Documents\MATLAB\github\RGC-Reconstruction\dat\wn_Dec3_sp_filter_win1.mat');


% [a,b]=which('rdtclient');
% 
% if isempty(a)
    load(filterFile);
    % if ismac || isunix
    %     load([reconstructionRootPath '/dat/' filterFile]);
    % else
    %     % load([phospheneRootPath '\dat\' filterFile]);
    %     load([reconstructionRootPath '\dat\' filterFile]);
    % end
    
% else
%     rdt = RdtClient('isetbio');
%     rdt.crp('/resources/data/rgc');
%     
%     data = rdt.readArtifact('pix1_nsBig_100Hz_1st_sv05__mosaicAll_8372855', 'type', 'mat');
%     filterMat = data.filterMat;
% end
% 
% load('C:\Users\James\Documents\MATLAB\github\RGC-Reconstruction\dat\wn_Dec3_sp_filter_new1.mat');

% rdt = RdtClient('isetbio');
% rdt.crp('/resources/data/reconstruction');
% data = rdt.readArtifact('filters_may26_parasol_midget_combined_svd_3000_len_100', 'type', 'mat');
% filterMat = data.filterMat; clear data;

% Reshape spikes and decode

recons_stim_on_off = cell(length(innerRetina.mosaic),1);
movrecons_on_off = cell(length(innerRetina.mosaic),1);

recons_stim_on_off_dropout = cell(length(innerRetina.mosaic),1);
movrecons_on_off_dropout = cell(length(innerRetina.mosaic),1);


blocklength = round(1/dt);
sizeSpikes = floor([ size(y{1},2) size(y{2},2) size(y{3},2) size(y{4},2)  ]/blocklength);

% on Parasol
spikesout = (y{1});%double(matfOff.spikesoutsm);
numcells1 = size(spikesout,1);

spikeRespOn = zeros(size(spikesout,1), max(sizeSpikes));

spikeRespOn(:,1:floor(size(spikesout,2)/blocklength)) = downSampResp(spikesout, numcells1, floor(size(spikesout,2)/blocklength),dt);
%     spikeRespOn= zeros(size(downSampResp(spikesout, numcells1, blocklength)));

% off Parasol
% numcells2 = 64;
spikesout2 = (y{2});%double(matfOff.spikesoutsm);
numcells2 = size(spikesout2,1);
spikeRespOff = zeros(size(spikesout2,1), max(sizeSpikes));
spikeRespOff(:,1:floor(size(spikesout2,2)/blocklength)) = ((downSampResp(spikesout2, numcells2, floor(size(spikesout2,2)/blocklength),dt)));
%     spikeRespOff = zeros(size(downSampResp(spikesout2, numcells2, blocklength)));

% off Midget
% numcells3 = 225;
spikesout3 = (y{3});%double(matfOff.spikesoutsm);

numcells3 = size(spikesout3,1);

spikeRespOffM = zeros(size(spikesout3,1), max(sizeSpikes));
spikeRespOffM(:,1:floor(size(spikesout3,2)/blocklength)) = ((downSampResp(spikesout3, numcells3, floor(size(spikesout3,2)/blocklength),dt)));
%     spikeRespOffM = zeros(size(downSampResp(spikesout3, numcells3, blocklength)));

% on Midget
% numcells4 = 169;
spikesout4 = (y{4});%double(matfOff.spikesoutsm);

numcells4 = size(spikesout4,1);
spikeRespOnM = zeros(size(spikesout4,1), max(sizeSpikes));
spikeRespOnM(:,1:floor(size(spikesout4,2)/blocklength)) = ((downSampResp(spikesout4, numcells4, floor(size(spikesout4,2)/blocklength),dt)));
%     spikeRespOnM = zeros(size(downSampResp(spikesout4, numcells4, blocklength)));

% spikesout = vertcat(onSR(:,1:15000), offSR(:,1:15000), onPSR(:,1:15000), offPSR(:,1:15000));

% whole mosaic
spikeRespOnOff =vertcat(spikeRespOn,spikeRespOff, spikeRespOffM,spikeRespOnM);
% only parasols
%     spikeRespOnOff =vertcat(spikeRespOn,spikeRespOn, zeros(size(spikeRespOffM)),zeros(size(spikeRespOnM)));
% only midgets
%     spikeRespOnOff =vertcat(zeros(size(spikeRespOn)),zeros(size(spikeRespOff)), ((spikeRespOffM)),((spikeRespOnM)));

% numbins = 4;
% spikeRespOnOff(1:890,:) = zeros(890,size(spikeRespOnOff,2));

% recons_stim_on_off{mosaicInd,1} = reconsFromFiltLen(filterMat, 1*spikeRespOnOff, numbins);

[mgr,mgc] = meshgrid(1:100,1:100);

[cmax,cind] = max(abs(filterMat),[],2);
[fmaxc,fmaxr] = ind2sub([100 100],cind);

mgrmat = mgr(:)*ones(1,size(fmaxr,1));
fmaxrmat = ones(size(mgrmat,1),1)*fmaxr';
mgrd = ((mgrmat - fmaxrmat)').^2;

mgcmat = mgc(:)*ones(1,size(fmaxc,1));
fmaxcmat = ones(size(mgcmat,1),1)*fmaxc';
mgcd = ((mgcmat - fmaxcmat)').^2;

dp = sqrt(mgrd+mgcd);
filterMat2 = filterMat;
filterMat2(dp>30) = 0;

recons_stim_on_off{mosaicInd,1} = reconsFromFiltLen(filterMat2, 1*spikeRespOnOff, numbins);

% % figure; ecdf(filterMat(:));
% filterMatInd = find(abs(filterMat)<0.0005); filterMat2 = filterMat; filterMat2(filterMatInd)=0;
% recons_stim_on_off{mosaicInd,1} = reconsFromFiltLen(filterMat2, 1*spikeRespOnOff, numbins);
szMov = sqrt(size(recons_stim_on_off{mosaicInd,1},1));
movrecons_on_off_full = reshape(recons_stim_on_off{mosaicInd,1},szMov,szMov,size(recons_stim_on_off{mosaicInd,1},2));

%% Reconstruction with a random subset of spiking cells zeroed out
deadIndicesAll = randperm(size(spikeRespOnOff,1));
numberDead = round(percentDead*(size(spikeRespOnOff,1)));
deadIndices = deadIndicesAll(1:numberDead);
spikeRespOnOffDead = spikeRespOnOff;
spikeRespOnOffDead(deadIndices,:) = zeros(length(deadIndices),size(spikeRespOnOff,2));

recons_stim_on_off_dropout{mosaicInd,1} = reconsFromFiltLen(filterMat, spikeRespOnOffDead, numbins);

movrecons_on_off_dropout = reshape(recons_stim_on_off_dropout{mosaicInd,1},szMov,szMov,size(recons_stim_on_off{mosaicInd,1},2));
% movrecons_on_off = .5*movrecons_on_off_full./mean(movrecons_on_off_full(:));
