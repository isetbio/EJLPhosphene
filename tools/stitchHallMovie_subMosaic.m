

parasolFlag = 0;
midgetFlag = 1;

load('/Users/james/Documents/MATLAB/EJLPhosphene/dat/stimuli/hallMovie.mat','vidFrame');

%%
load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/ws_hall_topleft_long.mat');
% movieStitch1 = movieStitch; movieStitchHealthy1 = movieStitchHealthy; electrodeArray1 = electrodeArray;
% clear movrecons_on_off movieStitch

if ~parasolFlag
    
    sp1 = innerRetinaHealthy.mosaic{1}.responseSpikes; %get('responseSpikes');
    size(sp1); sp1z = cell(size(sp1));
    for i = 1:size(sp1,1); for j = 1:size(sp1,2); sp1z{i,j} = []; end; end;
    innerRetinaHealthy.mosaic{1}.set('responseSpikes',sp1z);
    
    sp2 = innerRetinaHealthy2.mosaic{1}.get('responseSpikes');
    size(sp2); sp2z = cell(size(sp2));
    for i = 1:size(sp2,1); for j = 1:size(sp2,2); sp2z{i,j} = []; end; end;
    innerRetinaHealthy2.mosaic{1}.set('responseSpikes',sp2z);
end

if ~midgetFlag
    sp3 = innerRetinaHealthy3.mosaic{1}.responseSpikes; %get('responseSpikes');
    size(sp3); sp3z = cell(size(sp3));
    for i = 1:size(sp3,1); for j = 1:size(sp3,2); sp3z{i,j} = []; end; end;
    innerRetinaHealthy3.mosaic{1}.set('responseSpikes',sp3z);
    
    sp4 = innerRetinaHealthy4.mosaic{1}.get('responseSpikes');
    size(sp4); sp4z = cell(size(sp4));
    for i = 1:size(sp4,1); for j = 1:size(sp4,2); sp4z{i,j} = []; end; end;
    innerRetinaHealthy4.mosaic{1}.set('responseSpikes',sp4z);
end

[movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalRecon(innerRetinaHealthy, innerRetinaHealthy2, innerRetinaHealthy3, innerRetinaHealthy4, percentDead);
movieStitchHealthy1 = movrecons_on_offHealthy{1};

clear movieStitch movieStitchHealthy electrodeArray movrecons_on_offHealthy movrecons_on_offHealthy_dropout
%%
load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/ws_hall_bottomleft_long.mat');
% movieStitch2 = movieStitch; movieStitchHealthy2 = movieStitchHealthy; electrodeArray2 = electrodeArray;

if ~parasolFlag
    
    sp1 = innerRetinaHealthy.mosaic{1}.responseSpikes; %get('responseSpikes');
    size(sp1); sp1z = cell(size(sp1));
    for i = 1:size(sp1,1); for j = 1:size(sp1,2); sp1z{i,j} = []; end; end;
    innerRetinaHealthy.mosaic{1}.set('responseSpikes',sp1z);
    
    sp2 = innerRetinaHealthy2.mosaic{1}.get('responseSpikes');
    size(sp2); sp2z = cell(size(sp2));
    for i = 1:size(sp2,1); for j = 1:size(sp2,2); sp2z{i,j} = []; end; end;
    innerRetinaHealthy2.mosaic{1}.set('responseSpikes',sp2z);
end

if ~midgetFlag
    sp3 = innerRetinaHealthy3.mosaic{1}.responseSpikes; %get('responseSpikes');
    size(sp3); sp3z = cell(size(sp3));
    for i = 1:size(sp3,1); for j = 1:size(sp3,2); sp3z{i,j} = []; end; end;
    innerRetinaHealthy3.mosaic{1}.set('responseSpikes',sp3z);
    
    sp4 = innerRetinaHealthy4.mosaic{1}.get('responseSpikes');
    size(sp4); sp4z = cell(size(sp4));
    for i = 1:size(sp4,1); for j = 1:size(sp4,2); sp4z{i,j} = []; end; end;
    innerRetinaHealthy4.mosaic{1}.set('responseSpikes',sp4z);
end

[movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalRecon(innerRetinaHealthy, innerRetinaHealthy2, innerRetinaHealthy3, innerRetinaHealthy4, percentDead);
movieStitchHealthy2 = movrecons_on_offHealthy{1};

clear movieStitch movieStitchHealthy electrodeArray
%%
load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/ws_hall_topright_long.mat');
% movieStitch3 = movieStitch; movieStitchHealthy3 = movieStitchHealthy; electrodeArray3 = electrodeArray;

if ~parasolFlag
    
    sp1 = innerRetinaHealthy.mosaic{1}.responseSpikes; %get('responseSpikes');
    size(sp1); sp1z = cell(size(sp1));
    for i = 1:size(sp1,1); for j = 1:size(sp1,2); sp1z{i,j} = []; end; end;
    innerRetinaHealthy.mosaic{1}.set('responseSpikes',sp1z);
    
    sp2 = innerRetinaHealthy2.mosaic{1}.get('responseSpikes');
    size(sp2); sp2z = cell(size(sp2));
    for i = 1:size(sp2,1); for j = 1:size(sp2,2); sp2z{i,j} = []; end; end;
    innerRetinaHealthy2.mosaic{1}.set('responseSpikes',sp2z);
end

if ~midgetFlag
    sp3 = innerRetinaHealthy3.mosaic{1}.responseSpikes; %get('responseSpikes');
    size(sp3); sp3z = cell(size(sp3));
    for i = 1:size(sp3,1); for j = 1:size(sp3,2); sp3z{i,j} = []; end; end;
    innerRetinaHealthy3.mosaic{1}.set('responseSpikes',sp3z);
    
    sp4 = innerRetinaHealthy4.mosaic{1}.get('responseSpikes');
    size(sp4); sp4z = cell(size(sp4));
    for i = 1:size(sp4,1); for j = 1:size(sp4,2); sp4z{i,j} = []; end; end;
    innerRetinaHealthy4.mosaic{1}.set('responseSpikes',sp4z);
end


[movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalRecon(innerRetinaHealthy, innerRetinaHealthy2, innerRetinaHealthy3, innerRetinaHealthy4, percentDead);
movieStitchHealthy3 = movrecons_on_offHealthy{1};

clear movieStitch movieStitchHealthy electrodeArray
%%
load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/ws_hall_bottomright_long.mat');
% movieStitch4 = movieStitch; movieStitchHealthy4 = movieStitchHealthy; electrodeArray4 = electrodeArray;

if ~parasolFlag
    
    sp1 = innerRetinaHealthy.mosaic{1}.responseSpikes; %get('responseSpikes');
    size(sp1); sp1z = cell(size(sp1));
    for i = 1:size(sp1,1); for j = 1:size(sp1,2); sp1z{i,j} = []; end; end;
    innerRetinaHealthy.mosaic{1}.set('responseSpikes',sp1z);
    
    sp2 = innerRetinaHealthy2.mosaic{1}.get('responseSpikes');
    size(sp2); sp2z = cell(size(sp2));
    for i = 1:size(sp2,1); for j = 1:size(sp2,2); sp2z{i,j} = []; end; end;
    innerRetinaHealthy2.mosaic{1}.set('responseSpikes',sp2z);
end

if ~midgetFlag
    sp3 = innerRetinaHealthy3.mosaic{1}.responseSpikes; %get('responseSpikes');
    size(sp3); sp3z = cell(size(sp3));
    for i = 1:size(sp3,1); for j = 1:size(sp3,2); sp3z{i,j} = []; end; end;
    innerRetinaHealthy3.mosaic{1}.set('responseSpikes',sp3z);
    
    sp4 = innerRetinaHealthy4.mosaic{1}.get('responseSpikes');
    size(sp4); sp4z = cell(size(sp4));
    for i = 1:size(sp4,1); for j = 1:size(sp4,2); sp4z{i,j} = []; end; end;
    innerRetinaHealthy4.mosaic{1}.set('responseSpikes',sp4z);
end

[movrecons_on_offHealthy, movrecons_on_offHealthy_dropout] = irOptimalRecon(innerRetinaHealthy, innerRetinaHealthy2, innerRetinaHealthy3, innerRetinaHealthy4, percentDead);
movieStitchHealthy4 = movrecons_on_offHealthy{1};

clear movieStitch movieStitchHealthy electrodeArray


%%

movieStitchHealthy1rs = reshape(movieStitchHealthy1,[size(movieStitchHealthy1,1)*size(movieStitchHealthy1,2), size(movieStitchHealthy1,3)]);
m1 = mean(movieStitchHealthy1rs,1); 
movieStitchHealthy1rsmz = movieStitchHealthy1rs - repmat(mean(movieStitchHealthy1rs,1), [size(movieStitchHealthy1,1)*size(movieStitchHealthy1,2) 1]);
movieStitchHealthy1mz = reshape(movieStitchHealthy1rsmz,[size(movieStitchHealthy1,1),size(movieStitchHealthy1,2), size(movieStitchHealthy1,3)]);

movieStitchHealthy2rs = reshape(movieStitchHealthy2,[size(movieStitchHealthy2,1)*size(movieStitchHealthy2,2), size(movieStitchHealthy2,3)]);
m1 = mean(movieStitchHealthy2rs,1); 
movieStitchHealthy2rsmz = movieStitchHealthy2rs - repmat(mean(movieStitchHealthy2rs,1), [size(movieStitchHealthy2,1)*size(movieStitchHealthy2,2) 1]);
movieStitchHealthy2mz = reshape(movieStitchHealthy2rsmz,[size(movieStitchHealthy2,1),size(movieStitchHealthy2,2), size(movieStitchHealthy2,3)]);

movieStitchHealthy3rs = reshape(movieStitchHealthy3,[size(movieStitchHealthy3,1)*size(movieStitchHealthy3,2), size(movieStitchHealthy3,3)]);
m1 = mean(movieStitchHealthy3rs,1); 
movieStitchHealthy3rsmz = movieStitchHealthy3rs - repmat(mean(movieStitchHealthy3rs,1), [size(movieStitchHealthy3,1)*size(movieStitchHealthy3,2) 1]);
movieStitchHealthy3mz = reshape(movieStitchHealthy3rsmz,[size(movieStitchHealthy3,1),size(movieStitchHealthy3,2), size(movieStitchHealthy3,3)]);

movieStitchHealthy4rs = reshape(movieStitchHealthy4,[size(movieStitchHealthy4,1)*size(movieStitchHealthy4,2), size(movieStitchHealthy4,3)]);
m1 = mean(movieStitchHealthy4rs,1); 
movieStitchHealthy4rsmz = movieStitchHealthy4rs - repmat(mean(movieStitchHealthy4rs,1), [size(movieStitchHealthy4,1)*size(movieStitchHealthy4,2) 1]);
movieStitchHealthy4mz = reshape(movieStitchHealthy4rsmz,[size(movieStitchHealthy4,1),size(movieStitchHealthy4,2), size(movieStitchHealthy4,3)]);


%%
shiftval = 13;

% clear movieStitchHealthy
% movieStitchHealthy(1:96,1:96,:) = movieStitch1(:,:,1:500-shiftval) - mean(movieStitch1(:));
% movieStitchHealthy(96+[1:96],1:96,:) = movieStitch2(:,:,1:500-shiftval) - mean(movieStitch2(:));
% movieStitchHealthy(1:96,96+[1:96],:) = movieStitch3(:,:,1:500-shiftval) - mean(movieStitch3(:));
% movieStitchHealthy(96+[1:96],96+[1:96],:) = movieStitch4(:,:,1:500-shiftval) - mean(movieStitch4(:));

% % clear movieStitchHealthy
% healthyShift = 0;
% movieStitchHealthy(1+healthyShift:96,96+healthyShift+[1:96],:) = movieStitchHealthy1(:,:,1:500-shiftval) - mean(movieStitchHealthy1(:));
% movieStitchHealthy(96+healthyShift+[1:96],96+healthyShift+[1:96],:) = movieStitchHealthy2(:,:,1:500-shiftval) - mean(movieStitchHealthy2(:));
% movieStitchHealthy(1+healthyShift:96,96+healthyShift+96+[1:96],:) = movieStitchHealthy3(:,:,1:500-shiftval) - mean(movieStitchHealthy3(:));
% movieStitchHealthy(96+healthyShift+[1:96],96+healthyShift+96+[1:96],:) = movieStitchHealthy4(:,:,1:500-shiftval) - mean(movieStitchHealthy4(:));

% clear movieStitchHealthy
healthyShift = 0;
clear movieStitchHealthy
movieStitchHealthy(1+healthyShift:96,healthyShift+[1:96],:) = movieStitchHealthy1mz(:,:,1:500-shiftval);% - mean(movieStitchHealthy1(:));
movieStitchHealthy(96+healthyShift+[1:96],healthyShift+[1:96],:) = movieStitchHealthy2mz(:,:,1:500-shiftval);% - mean(movieStitchHealthy2(:));
movieStitchHealthy(1+healthyShift:96,healthyShift+96+[1:96],:) = movieStitchHealthy3mz(:,:,1:500-shiftval);% - mean(movieStitchHealthy3(:));
movieStitchHealthy(96+healthyShift+[1:96],healthyShift+96+[1:96],:) = movieStitchHealthy4mz(:,:,1:500-shiftval);% - mean(movieStitchHealthy4(:));

vf = vidFrame(:,:,shiftval:500);
stimShift = 0;
movieStitchHealthy([1:192],stimShift+192+[1:192],1:500-shiftval+1) = (vf-mean(vf(:)))./25;

% load('/Users/james/Documents/MATLAB/EJLPhosphene/dat/movieStitchHealthyFull.mat','movieStitchHealthyFull');
load('/Users/james/Documents/MATLAB/EJLPhosphene/dat/movieStitchHealthyFull.mat','movieStitchHealthyFull');
movieStitchHealthy = movieStitchHealthyFull;
figure; 
set(gcf,'position',[483         703        1077         635]);
ieMovie(movieStitchHealthy);
% figure; ieMovie(movieStitchHealthy,'save',true,'vname','protheseis_navigation2');

%%
% figure; 
% subplot(121);
% imagesc(movieStitchHealthy(:,:,450));
% subplot(122);
% imagesc(vidFrame(:,:,486));

%%

