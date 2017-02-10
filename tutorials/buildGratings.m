function gratingsMovie = buildGratings(numRows,numCols,gratingFreq)
% Assume sampling rate is 30 Hz

% intervalSteps = 10;
% telec = 2;

gratingsMovie = 0.5*ones(numRows,numCols);
telec = 0;
% gratingsMovie(:) = 0;

% % moving bar
% for telec = 1:25
%     gratingsMovie(telec,:,5*(telec-1)+1:5*telec) = 1;
% end

% full field
intervalSteps = gratingFreq;

telectot = 0;
for gratingSpFreq =  [2.^(1:5) 50]
    for telec = 1:8*4
%     telec = telec+2;
    gratingsMovieRow = sin((2*pi/(numCols/gratingSpFreq))*[1:numCols]);
    tempFreq = sin(2*pi*telec/(8));
    
    gratingsMovie(:,:,telectot+telec) = ones(numRows,1)*gratingsMovieRow*tempFreq;
    end
    telectot = telectot+telec;
end

% figure; ieMovie(gratingsMovie);

% for gratingSpFreq =  [2.^(1:7)]
% for gratingsRep = 1:4;
%     
%     telec = telec+2;
%     %     Full field
%     %     gratingsMovie(:,:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% 
%     % Half Field
%     gratingsMovieRow = sin((2*pi/(numCols/gratingSpFreq))*[1:numCols]);
%     for colStart = 1:2*round(numCols/gratingSpFreq):numCols
% %         gratingsMovie(:,1:round(numRows/gratingSpFreq)-1,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
%         gratingsMovie(:,colStart+[1:round(numRows/gratingSpFreq)-1]-1,intervalSteps*(telec-1)+1:intervalSteps*telec) = -.5*((-1)^(gratingsRep))+.5;
%         gratingsMovie(:,round(numCols/gratingSpFreq)+colStart+[1:round(numRows/gratingSpFreq)-1]-1,intervalSteps*(telec-1)+1:intervalSteps*telec) = -.5*((-1)^(gratingsRep+1))+.5;
%         
%     end
% end
% end
% figure; ieMovie(gratingsMovie);