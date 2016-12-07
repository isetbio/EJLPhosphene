function [movrecons_on_off_full, movrecons_on_off_dropout] = irOptimalRecon(innerRetina, innerRetina2, innerRetina3, innerRetina4, percentDead)

%%
y = cell(length(innerRetina.mosaic));
for mosaicInd = 1:length(innerRetina.mosaic)
    
    cellCtr=0; dt = .01;
    maxTrials = innerRetina.mosaic{1}.numberTrials;
    nCells = size(innerRetina.mosaic{1}.responseSpikes);
    %         yout = [];
    y{mosaicInd} = zeros(36,10000);
    for ycell = 1:nCells(1)
        for xcell = 1:nCells(2)
            %         clear yind y
            cellCtr = cellCtr+1;
            
            for trial = 1:maxTrials
                
                yind =  innerRetina.mosaic{mosaicInd}.responseSpikes{xcell,ycell,trial,1};
                
                %             y(xcell,ycell,trial,ceil(yind./dt))=1;
                y{mosaicInd}(cellCtr,ceil(yind./dt))=1;
                
            end
        end
    end
    
end

%%

y2 = cell(length(innerRetina2.mosaic));
for mosaicInd = 1:length(innerRetina2.mosaic)

cellCtr=0; dt = .01;
maxTrials = innerRetina2.mosaic{1}.numberTrials;
nCells = size(innerRetina2.mosaic{1}.responseSpikes);
%         yout = [];
y2{mosaicInd} = zeros(64,10000);
for ycell = 1:nCells(1)
    for xcell = 1:nCells(2)
%         clear yind y
        cellCtr = cellCtr+1;
        
        for trial = 1:maxTrials
            
            yind2 =  innerRetina2.mosaic{mosaicInd}.responseSpikes{xcell,ycell,trial,1};

%             y(xcell,ycell,trial,ceil(yind./dt))=1;
            y2{mosaicInd}(cellCtr,ceil(yind2./dt))=1;

        end
    end
end

end

%%


y3 = cell(length(innerRetina3.mosaic));
for mosaicInd = 1:length(innerRetina3.mosaic)

cellCtr=0; dt = .01;
maxTrials = innerRetina3.mosaic{1}.numberTrials;
nCells = size(innerRetina3.mosaic{1}.responseSpikes);
%         yout = [];
y3{mosaicInd} = zeros(225,10000);
for ycell = 1:nCells(1)
    for xcell = 1:nCells(2)
%         clear yind y
        cellCtr = cellCtr+1;
        
        for trial = 1:maxTrials
            
            yind3 =  innerRetina3.mosaic{mosaicInd}.responseSpikes{xcell,ycell,trial,1};

%             y(xcell,ycell,trial,ceil(yind./dt))=1;
            y3{mosaicInd}(cellCtr,ceil(yind3./dt))=1;

        end
    end
end

end

%%


y4 = cell(length(innerRetina4.mosaic));
for mosaicInd = 1:length(innerRetina4.mosaic)

cellCtr=0; dt = .01;
maxTrials = innerRetina4.mosaic{1}.numberTrials;
nCells = size(innerRetina4.mosaic{1}.responseSpikes);
%         yout = [];
y4{mosaicInd} = zeros(144,10000);
for ycell = 1:nCells(1)
    for xcell = 1:nCells(2)
%         clear yind y
        cellCtr = cellCtr+1;
        
        for trial = 1:maxTrials
            
            yind4 =  innerRetina4.mosaic{mosaicInd}.responseSpikes{xcell,ycell,trial,1};

%             y(xcell,ycell,trial,ceil(yind./dt))=1;
            y4{mosaicInd}(cellCtr,ceil(yind4./dt))=1;

        end
    end
end

end

%%

% sep 25
% load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/filters_may26_parasol_midget_combined_svd_3000_len_100.mat');
load('C:\Users\James\Documents\GitHub\RGC-ReconAdd\output\svd_reconstruct_shorttrain_midgets\filters_may26_parasol_midget_combined_svd_3000_len_100.mat')

% load('C:\Users\James\Documents\MATLAB\github\RGC-Reconstruction\dat\ns_Dec2_sp_filter.mat');
% rdt = RdtClient('isetbio');
% rdt.crp('/resources/data/reconstruction');
% data = rdt.readArtifact('filters_may26_parasol_midget_combined_svd_3000_len_100', 'type', 'mat');
% filterMat = data.filterMat; clear data;

% oct 8
% load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/oct8/filters_may26_fine_all2_svd_3000_len_100.mat')
% load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/oct8/filters_may26_fine_all2_svd_1500_len_100.mat')

% load('C:\Users\James\Documents\GitHub\RGC-ReconAdd\output\svd_reconstruct_shorttrain_midgets\filters_may26_on_long300_shorttime_svd_1000_len_100.mat')
%     numbins = 12;
%     recons_stim_on_off = reconsFromFiltLen(filterMat, spikeRespOn, numbins);
    
%     mov = reshape(stim,96,96,size(stim,2));
% %     movrecons_on_off = reshape(recons_stim_on_off,96,96,size(recons_stim_on_off,2));
%         movrecons_on_off_full = reshape(recons_stim_on_off,96,96,size(recons_stim_on_off,2));
%     movrecons_on_off =(.25+.5*(movrecons_on_off_full)./median((movrecons_on_off_full(:))));
%%




%
recons_stim_on_off = cell(length(innerRetina2.mosaic),1);
movrecons_on_off = cell(length(innerRetina2.mosaic),1);

recons_stim_on_off_dropout = cell(length(innerRetina2.mosaic),1);
movrecons_on_off_dropout = cell(length(innerRetina2.mosaic),1);
for mosaicInd = 1:length(innerRetina.mosaic)
    
    blocklength = 100;
    sizeSpikes = floor([ size(y{mosaicInd},2) size(y2{mosaicInd},2) size(y3{mosaicInd},2) size(y4{mosaicInd},2)  ]/blocklength);
    
    % on Parasol
    spikesout = (y{mosaicInd});%double(matfOff.spikesoutsm);
    numcells1 = 36;
%     downSampRespSetLength

    spikeRespOn = zeros(size(spikesout,1), max(sizeSpikes));

    spikeRespOn(:,1:floor(size(spikesout,2)/blocklength)) = downSampResp(spikesout, numcells1, floor(size(spikesout,2)/blocklength));
%     spikeRespOn= zeros(size(downSampResp(spikesout, numcells1, blocklength)));
    
    % off Parasol
    numcells2 = 64;
    spikesout2 = (y2{mosaicInd});%double(matfOff.spikesoutsm);
    spikeRespOff = zeros(size(spikesout2,1), max(sizeSpikes));
    spikeRespOff(:,1:floor(size(spikesout2,2)/blocklength)) = ((downSampResp(spikesout2, numcells2, floor(size(spikesout2,2)/blocklength))));
%     spikeRespOff = zeros(size(downSampResp(spikesout2, numcells2, blocklength)));

    % off Midget
    numcells3 = 225;
    spikesout3 = (y3{mosaicInd});%double(matfOff.spikesoutsm);
    
    spikeRespOffM = zeros(size(spikesout3,1), max(sizeSpikes));
    spikeRespOffM(:,1:floor(size(spikesout3,2)/blocklength)) = ((downSampResp(spikesout3, numcells3, floor(size(spikesout3,2)/blocklength))));
%     spikeRespOffM = zeros(size(downSampResp(spikesout3, numcells3, blocklength)));

    % on Midget
    numcells4 = 144;
    spikesout4 = (y4{mosaicInd});%double(matfOff.spikesoutsm);
    
    spikeRespOnM = zeros(size(spikesout4,1), max(sizeSpikes));
    spikeRespOnM(:,1:floor(size(spikesout4,2)/blocklength)) = ((downSampResp(spikesout4, numcells4, floor(size(spikesout4,2)/blocklength))));
%     spikeRespOnM = zeros(size(downSampResp(spikesout4, numcells4, blocklength)));
        
    % spikesout = vertcat(onSR(:,1:15000), offSR(:,1:15000), onPSR(:,1:15000), offPSR(:,1:15000));
    
    % whole mosaic
    spikeRespOnOff =vertcat(spikeRespOnM,spikeRespOffM, spikeRespOn,spikeRespOff);
    % only midgets
%     spikeRespOnOff =vertcat(spikeRespOnM,spikeRespOffM, zeros(size(spikeRespOn)),zeros(size(spikeRespOff)));
    % only parasols
%     spikeRespOnOff =vertcat(zeros(size(spikeRespOnM)),zeros(size(spikeRespOffM)), ((spikeRespOn)),((spikeRespOff)));
    
    numbins = 8;
      filterMatInd = find(abs(filterMat)<0.002); filterMat2 = filterMat; filterMat2(filterMatInd)=0;
    recons_stim_on_off{mosaicInd,1} = reconsFromFiltLen(filterMat2, spikeRespOnOff, numbins);      
    movrecons_on_off_full{mosaicInd,1} = reshape(recons_stim_on_off{mosaicInd,1},96,96,size(recons_stim_on_off{mosaicInd,1},2));  
    
    deadIndicesAll = randperm(size(spikeRespOnOff,1));
    numberDead = round(percentDead*(size(spikeRespOnOff,1)));
    deadIndices = deadIndicesAll(1:numberDead);
    spikeRespOnOffDead = spikeRespOnOff;
    spikeRespOnOffDead(deadIndices,:) = zeros(length(deadIndices),size(spikeRespOnOff,2));
        
    recons_stim_on_off_dropout{mosaicInd,1} = reconsFromFiltLen(filterMat, spikeRespOnOffDead, numbins);
      
    movrecons_on_off_dropout{mosaicInd,1} = reshape(recons_stim_on_off_dropout{mosaicInd,1},96,96,size(recons_stim_on_off{mosaicInd,1},2));
    % movrecons_on_off = .5*movrecons_on_off_full./mean(movrecons_on_off_full(:));
    
end