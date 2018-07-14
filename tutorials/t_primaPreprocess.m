% t_primaPreprocess
% 
% The linear model for the healthy retina is
% 
%       im = W*sp           (decoding)
%       sp = W^(-1)*im      (encoding)
% 
% The linear model for the prothesis stimulation is
% 
%       im_p = W_p*sp_p         (decoding)
%       sp_p = W_p^(-1)*im_p    (encoding)
% 
% If we have the mismatch case where
% 
%       im_mis = W*sp_p = W*(W_p^(-1)*im_p)
%       im_p = W^(-1)*Wp*im_mis
% 
% Thus the ideal linear preprocessing transform is W^(-1)*Wp.
% 
% This does not work yet!

% 
%% Load the healthy decoding filter and take pinv(W)

clear
% load('filters_mosaic0_sv80_w1_sh4_may22.mat')
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/current/filters_mosaic0_sv20_w1_sh2_dr0.mat')
 
load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv50_w1_sh15_dr0.mat');
filterNS = filterMat; clear filterMat;
% rd = RdtClient('isetbio');
% rd.crp('/resources/data/istim');
% filterFile = 'filters_mosaic0_sv80_w1_sh4_may22.mat';
% data  = rd.readArtifact(filterFile(1:end-4), 'type', 'mat');
% filterMat = data.filterMat; clear data;

% Apply the exponential zeroing filter W_smoothed
filterHealthyExp = zeroFilter(filterNS,.06);
% filterHealthyExp = zeroFilter(filterMat,.0075);
figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(filterHealthyExp(2400+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;

% filterHealthyExpsc = filterHealthyExp./((ones(size(filterHealthyExp,2),1)*max(abs(filterHealthyExp),[],2)')');
filterHealthyExpsc = filterHealthyExp./ (sqrt(ones(size(filterHealthyExp,2),1)*(sum(filterHealthyExp.^2,2))')');


% filterHealthyExpzm = (filterHealthyExp-(ones(size(filterHealthyExp,2),1)*mean(filterHealthyExp,2)')');
% filterHealthyExpsc = filterHealthyExpzm./((ones(size(filterHealthyExpzm,2),1)*max(abs(filterHealthyExpzm),[],2)')');
% 
% % Visualize
% figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(filterHealthyExpsc(2400+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;
% % absFilter = abs(filterHealthyExp);


%%

% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv50_w1_sh15_dr0.mat');
% filterNS = filterMat; clear filterMat
% load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv 5_w1_sh4_dr0_pitch_70_decay_2.mat')
% filterPros = filterMat;

%%
  load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv 5_w1_sh4_dr0_pitch_70_decay_2.mat')
  filterPros = filterMat; clear filterMat;
    filterMatZero = zeroFilter(filterPros,.005);
    filterMatZerosc = filterMatZero./ (sqrt(ones(size(filterMatZero,2),1)*(sum(filterMatZero.^2,2))')');

    wopt = filterHealthyExpsc\filterMatZerosc;
    wopt2 = filterMatZerosc\filterHealthyExpsc;
    
    woptz = zeroFilter(wopt,.0005);
%

%%

load('/Users/james/Documents/MATLAB/RGC-Reconstruction/dat/aug29prima70/filtersmosaic0_sv 2_w1_sh3_dr0.mat')
filter35 = filterMat; clear filterMat;
   filter35Zero = zeroFilter(filter35,.005);
    filter35Zerosc = filter35Zero./ (sqrt(ones(size(filter35Zero,2),1)*(sum(filter35Zero.^2,2))')');

    wopt = filterHealthyExpsc\filter35Zerosc;
  

%% Take pinv(W_smoothed)
filterHealthyExpsc = zeroFilter(filterNS,.01);
pinvFilterHealthyExp = pinv(filterHealthyExpsc);

figure; imagesc(reshape(pinvFilterHealthyExp(:,100),[100 100]))

save('pinvHealthyZero2.mat','pinvFilterHealthyExp');

%% Load the prosthesis decoding filter and tkae pinv(Wp)

% load('filters_mosaic0_sv75_w1_sh2_may25prima.mat')
% load('filters_mosaic0_sv75_w1_sh2_may26primaSmall.mat')

rd = RdtClient('isetbio');
rd.crp('/resources/data/istim');
% filterFile = 'filters_mosaic0_sv75_w1_sh2_may26primaSmall';
filterFile = 'filters_mosaic0_sv20_w1_sh2_dr0';
data  = rd.readArtifact(filterFile, 'type', 'mat');
filterMat = data.filterMat; clear data;

% Smoothing filter to generate Wp_smoothed
filterMatZero = zeroFilter(filterPros,.001);
filterMatZerozm = (filterMatZero-(ones(size(filterMatZero,2),1)*mean(filterMatZero,2)')');
filterMatZerosc = filterMatZerozm./((ones(size(filterMatZerozm,2),1)*max(abs(filterMatZerozm),[],2)')');

figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(filterMat(100+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;

figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(filterMatZerosc(100+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;


figure; imagesc(reshape(filterMatZerosc(100,:),[100 100]))
% 
% noLearningImageFilter = pinvFilterHealthy*filterMatZero;

noLearningImageFilter = pinvFilterHealthyExp*filterMatZero;

noLearningImageFilter2 = filterMatZero'*pinvFilterHealthyExp';
% 
% % figure; imagesc(noLearningImageFilter(1:500,:))
% % figure; imagesc(reshape(noLearningImageFilter(:,100),[100 100]))
% % figure; imagesc(reshape(noLearningImageFilter(100,:),[100 100]))
% % figure; imagesc(reshape(noLearningImageFilter(1001,:),[100 100]))
% % figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(noLearningImageFilter(2400+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;
% % figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(noLearningImageFilter(100+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;
% % figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(noLearningImageFilter(700+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;
% % figure; imagesc(reshape(sum(noLearningImageFilter(:,:).^2),[100 100]))
% 
% 
% 
%% Load stimulus movie

stimSize = 100;

rsFactor = 1;

load([phospheneRootPath '/dat/stimuli/hallMovie.mat']);
szFrames = size(vidFrame,3);
hallMovieResize = zeros(rsFactor*stimSize,rsFactor*stimSize,szFrames);
for ii = 1:szFrames
    hallMovieResize(:,:,ii) = imresize(vidFrame(:,:,ii),[rsFactor*stimSize,rsFactor*stimSize]);
end

% movieIn = rand(100,100,10);
movieIn = hallMovieResize(:,:,1:600);

clear vidFrame hallMovieResize

movieFr = movieIn(:,:,100);

movieFrFilt = noLearningImageFilter*movieFr(:);
figure; imagesc(reshape(movieFrFilt,[100 100]))

%% Try applying the filter and doing some zeroing

noLearningImageFilterzm = (noLearningImageFilter-(ones(size(noLearningImageFilter,2),1)*mean(noLearningImageFilter,2)')');
noLearningImageFiltersc = noLearningImageFilterzm./((ones(size(noLearningImageFilterzm,2),1)*max(noLearningImageFilterzm,[],2)')');

figure; imagesc(noLearningImageFiltersc(1:100,:))
figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(noLearningImageFiltersc(700+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;


movieFrFilt = noLearningImageFiltersc'*movieFr(:);
figure; imagesc(reshape(movieFrFilt,[100 100]))