


clear
load('filters_mosaic0_sv80_w1_sh4_may22.mat')
 
filterHealthyExp = zeroFilter(filterMat,.0075);
filterHealthyExpzm = (filterHealthyExp-(ones(size(filterHealthyExp,2),1)*mean(filterHealthyExp,2)')');
filterHealthyExpsc = filterHealthyExpzm./((ones(size(filterHealthyExpzm,2),1)*max(abs(filterHealthyExpzm),[],2)')');


figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(filterHealthyExpsc(2400+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;
% absFilter = abs(filterHealthyExp);

%%
pinvFilterHealthyExp = pinv(filterHealthyExpsc);

figure; imagesc(reshape(pinvFilterHealthyExp(:,100),[100 100]))

save('pinvHealthyZero2.mat','pinvFilterHealthyExp');

%%
load('filters_mosaic0_sv75_w1_sh2_may25prima.mat')
% load('filters_mosaic0_sv75_w1_sh2_may26primaSmall.mat')

filterMatZero = zeroFilter(filterMat,.01);
filterMatZerozm = (filterMatZero-(ones(size(filterMatZero,2),1)*mean(filterMatZero,2)')');
filterMatZerosc = filterMatZerozm./((ones(size(filterMatZerozm,2),1)*max(abs(filterMatZerozm),[],2)')');

figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(filterMat(100+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;

figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(filterMatZerosc(100+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;


figure; imagesc(reshape(filterMatZerosc(100,:),[100 100]))
% 
% noLearningImageFilter = pinvFilterHealthy*filterMatZero;

noLearningImageFilter = pinvFilterHealthyExp*filterMatZero;

noLearningImageFilter2 = filterMatZerosc'*pinvFilterHealthyExp';
% 
% % figure; imagesc(noLearningImageFilter(1:500,:))
% % figure; imagesc(reshape(noLearningImageFilter(:,100),[100 100]))
% % figure; imagesc(reshape(noLearningImageFilter(100,:),[100 100]))
% % figure; imagesc(reshape(noLearningImageFilter(1001,:),[100 100]))
% % figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(noLearningImageFilter(2400+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;
% % figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(noLearningImageFilter(400+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;
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

movieFrFilt = noLearningImageFilter(:,100)'*movieFr(:);
figure; imagesc(reshape(movieFrFilt,[100 100]))

%%

noLearningImageFilterzm = (noLearningImageFilter-(ones(size(noLearningImageFilter,2),1)*mean(noLearningImageFilter,2)')');
noLearningImageFiltersc = noLearningImageFilterzm./((ones(size(noLearningImageFilterzm,2),1)*max(noLearningImageFilterzm,[],2)')');

figure; imagesc(noLearningImageFiltersc(1:100,:))
figure; for fr = 1:25;     subplot(5,5,fr);     filtIm = reshape(noLearningImageFiltersc(700+fr,:),[100 100]);     imagesc(filtIm);     mabs = max(abs(filtIm(:)));     caxis([-mabs mabs]); colormap parula;  end;


movieFrFilt = noLearningImageFiltersc'*movieFr(:);
figure; imagesc(reshape(movieFrFilt,[100 100]))