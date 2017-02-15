function innerRetina = irGetLinearRespElectrode(innerRetina, innerRetinaInput, innerRetinaThreshold, params)



for mosaicInd = 1:length(innerRetina.mosaic)
    clear innerRetinaActivation  i0all xc yc
    [yc xc] = size(innerRetina.mosaic{mosaicInd}.cellLocation);
    i0all{mosaicInd} = .5*median(innerRetinaInput(:))*ones(xc,yc);% -.05 - .45*rand(xc,yc);
%     tc = innerRetina.mosaic{1}.tCenter{1};
    for xind = 1:xc
        for yind = 1:yc
            
%             for frame = 1:params.nSteps
                % innerRetinaActivation{xind,yind,mosaicInd} = innerRetinaFunction{xind,yind,mosaicInd}(innerRetinaInput(xind,yind,mosaicInd));
                funcParams = innerRetinaThreshold{xind,yind,mosaicInd};
                thr = 80;%funcParams(1); 
                i0 = i0all{mosaicInd}(xind,yind);
                % innerRetinaFunction = @(iElectrode) 10./(1+exp(-thr*(i0+iElectrode)));
                % innerRetinaFunction = @(iElectrode) log(1./(1+exp(-thr*(i0+iElectrode))));
                
                % innerRetinaActivation{xind,yind,mosaicInd} = innerRetinaFunction(innerRetinaInput(xind,yind,mosaicInd));
                
%                 innerRetinaFunction =  @(iElectrode) (50*(1./(1+exp(-thr*(i0+iElectrode)))));
%                 innerRetinaActivation{xind,yind}(frame) = innerRetinaFunction(innerRetinaInput(xind,yind,frame,mosaicInd));
                
                innerRetinaFunction = @(iElectrode) (5*iElectrode);
                % innerRetinaActivation{xind,yind}(frame) = 50*innerRetinaFunction(innerRetinaInput(xind,yind,frame,mosaicInd));
                innerRetinaActivation(xind,yind,:) = 50*innerRetinaFunction(innerRetinaInput(xind,yind,:,mosaicInd));
%                 innerRetinaActivation(xind,yind,:) = conv(squeeze(50*innerRetinaFunction(innerRetinaInput(xind,yind,:,mosaicInd))),tc,'full');
%             end
        end
    end
    innerRetinaActivation0{mosaicInd} = innerRetinaActivation;
    
    innerRetinaActivationZM = innerRetinaActivation-mean(innerRetinaActivation(:));
    
    % Threshold using ECDF - knock lowest 10% of abs(activations) to zero
    [fe,xe]=ecdf(abs(innerRetinaActivationZM(:))); % figure; plot(xe,fe)
    threshECDF = 0.1; % set threshold level
    threshVal = xe(min(find(fe>threshECDF)));    
    
    mosaicSet(innerRetina.mosaic{mosaicInd},'responseLinear', (innerRetinaActivationZM.*innerRetinaActivationZM>threshVal));
    % mosaicSet(innerRetina2.mosaic{mosaicInd},'responseLinear', innerRetinaActivation);

end

% irPlot(innerRetina, 'linear');