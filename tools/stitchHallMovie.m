


load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/hallMovie.mat','vidFrame');

load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/ws_hall_topleft_long.mat');
movieStitch1 = movieStitch; movieStitchHealthy1 = movieStitchHealthy; electrodeArray1 = electrodeArray;
clear movieStitch movieStitchHealthy electrodeArray
load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/ws_hall_bottomleft_long.mat');
movieStitch2 = movieStitch; movieStitchHealthy2 = movieStitchHealthy; electrodeArray2 = electrodeArray;
clear movieStitch movieStitchHealthy electrodeArray
load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/ws_hall_topright_long.mat');
movieStitch3 = movieStitch; movieStitchHealthy3 = movieStitchHealthy; electrodeArray3 = electrodeArray;
clear movieStitch movieStitchHealthy electrodeArray
load('/Users/james/Documents/MATLAB/isetbio misc/eye_and_chip/sep25/ws_hall_bottomright_long.mat');
movieStitch4 = movieStitch; movieStitchHealthy4 = movieStitchHealthy; electrodeArray4 = electrodeArray;
clear movieStitch movieStitchHealthy electrodeArray


%%
shiftval = 13;

clear movieStitchHealthy
movieStitchHealthy(1:96,1:96,:) = movieStitch1(:,:,1:500-shiftval) - mean(movieStitch1(:));
movieStitchHealthy(96+[1:96],1:96,:) = movieStitch2(:,:,1:500-shiftval) - mean(movieStitch2(:));
movieStitchHealthy(1:96,96+[1:96],:) = movieStitch3(:,:,1:500-shiftval) - mean(movieStitch3(:));
movieStitchHealthy(96+[1:96],96+[1:96],:) = movieStitch4(:,:,1:500-shiftval) - mean(movieStitch4(:));

% clear movieStitchHealthy
movieStitchHealthy(1:96,192+[1:96],:) = movieStitchHealthy1(:,:,1:500-shiftval) - mean(movieStitchHealthy1(:));
movieStitchHealthy(96+[1:96],192+[1:96],:) = movieStitchHealthy2(:,:,1:500-shiftval) - mean(movieStitchHealthy2(:));
movieStitchHealthy(1:96,192+96+[1:96],:) = movieStitchHealthy3(:,:,1:500-shiftval) - mean(movieStitchHealthy3(:));
movieStitchHealthy(96+[1:96],192+96+[1:96],:) = movieStitchHealthy4(:,:,1:500-shiftval) - mean(movieStitchHealthy4(:));

vf = vidFrame(:,:,shiftval:500);
movieStitchHealthy([1:192],192+192+[1:192],1:500-shiftval+1) = (vf-mean(vf(:)))./25;


figure; ieMovie(movieStitchHealthy);
% figure; ieMovie(movieStitchHealthy,'save',true,'vname','protheseis_navigation2');

%%
% figure; 
% subplot(121);
% imagesc(movieStitchHealthy(:,:,450));
% subplot(122);
% imagesc(vidFrame(:,:,486));

%%

