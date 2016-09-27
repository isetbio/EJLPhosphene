function innerRetinaThreshold = irGetThreshold(innerRetina)

for mosaicInd = 1:length(innerRetina.mosaic)
    [yc xc] = size(innerRetina.mosaic{mosaicInd}.cellLocation);
    for xind = 1:xc
        for yind = 1:yc
            thr = 20*rand(1,1);
            i0 = -.5*rand(1,1);
            % innerRetinaFunction{xind,yind,mosaicInd}= @(iElectrode) 1./(1+exp(-thr*(i0+iElectrode)));
            innerRetinaThreshold{xind,yind,mosaicInd} = [thr i0];
            % plot( innerRetinaFunction{xind,yind,mosaicInd}(-1:.1:4));
        end
    end
end