function fullFieldMovie = buildFullField(numRows,numCols)


intervalSteps = 10;
telec = 2;

fullFieldMovie = zeros(numRows,numCols);
% fullFieldMovie(:) = 0;

% % moving bar
% for telec = 1:25
%     fullFieldMovie(telec,:,5*(telec-1)+1:5*telec) = 1;
% end

% full field
intervalSteps = 10;
telec = 2;
fullFieldMovie(:,:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

% half field 1
telec = telec+2;
fullFieldMovie(:,1:round(numRows/2)-1,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

% % half field 2
telec = telec+2;
fullFieldMovie(:,round(numRows/2):end,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

% % half field 3
% telec = telec+2;
% fullFieldMovie(1:12,:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

% % half field 4
% telec = telec+2;
% fullFieldMovie(13:end,:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

%  quarter field 1
telec = telec+2;
fullFieldMovie(round(numCols/2):end,1:round(numRows/2)-1,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

%  quarter field 2
telec = telec+2;
fullFieldMovie(1:round(numCols/2)-1,round(numRows/2):end,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

xex = ceil(numRows/30); yex = ceil(numCols/30);

% subsets of six at four points
% point 1
telec = telec+2;
centerv = [round(numCols/3) round(numRows/3)]; hexcenter = [];
cctr = 0; for c1 = -xex:xex; for c2 = -yex:yex; cctr = cctr+1; hexcenter(cctr,:) = [centerv(1)+c1 centerv(2)+c2]; end; end;
fullFieldMovie(hexcenter(:,1),hexcenter(:,2),intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% fullFieldMovie([1:2 end-1:end],:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% fullFieldMovie(:,[1:2 end-1:end],intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% point 2
centerv = [round(numCols/3) round(2*numRows/3)]; hexcenter = [];
cctr = 0; for c1 = -xex:xex; for c2 = -yex:yex; cctr = cctr+1; hexcenter(cctr,:) = [centerv(1)+c1 centerv(2)+c2]; end; end;
fullFieldMovie(hexcenter(:,1),hexcenter(:,2),intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% fullFieldMovie([1:2 end-1:end],:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% fullFieldMovie(:,[1:2 end-1:end],intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% point 3
centerv = [round(2*numCols/3)-1 round(numRows/3)]; hexcenter = [];
cctr = 0; for c1 = -xex:xex; for c2 = -yex:yex; cctr = cctr+1; hexcenter(cctr,:) = [centerv(1)+c1 centerv(2)+c2]; end; end;
fullFieldMovie(hexcenter(:,1),hexcenter(:,2),intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% fullFieldMovie([1:2 end-1:end],:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% fullFieldMovie(:,[1:2 end-1:end],intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% point 4
centerv = [round(2*numCols/3)-1 round(2*numRows/3)]; hexcenter = [];
cctr = 0; for c1 = -xex:xex; for c2 = -yex:yex; cctr = cctr+1; hexcenter(cctr,:) = [centerv(1)+c1 centerv(2)+c2]; end; end;
fullFieldMovie(hexcenter(:,1),hexcenter(:,2),intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% fullFieldMovie([1:2 end-1:end],:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;
% fullFieldMovie(:,[1:2 end-1:end],intervalSteps*(telec-1)+1:intervalSteps*telec) = 1;

telec = telec+2;
fullFieldMovie(:,:,intervalSteps*(telec-1)+1:intervalSteps*telec) = 0;


%%%% Border activation
xborder = 4; yborder = 3;
fullFieldMovie([1:xborder end-xborder+1:end],:,:) = .25;
fullFieldMovie(:,[1:yborder end-yborder+1:end],:) = .25;
