function innerRetina = scaleRF(innerRetina, scaleFactor)

for mosaicInd = 1:length(innerRetina.mosaic)
    for ri = 1:size((innerRetina.mosaic{mosaicInd}.cellLocation),1);
        for ci = 1:size((innerRetina.mosaic{mosaicInd}.cellLocation),2);
            sRFcenterNew{ri,ci} = scaleFactor*innerRetina.mosaic{mosaicInd}.sRFcenter{ri,ci};
            sRFsurroundNew{ri,ci} = scaleFactor*innerRetina.mosaic{mosaicInd}.sRFsurround{ri,ci};
        end
    end
    innerRetina.mosaic{mosaicInd}.set('sRFcenter', sRFcenterNew);
    innerRetina.mosaic{mosaicInd}.set('sRFsurround', sRFsurroundNew);
end