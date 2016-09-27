

v = VideoReader('/Users/james/Downloads/IMG_4056.MOV');

       % Specify that reading should start at 0.5 seconds from the
        % beginning.
        v.CurrentTime = 0.5;
 
        % Create an axes
        currAxes = axes;
        
        % Read video frames until available
        frctr = 0;
        vidFrame = zeros(192,192,30*20);
        while hasFrame(v)
            clear frameTemp
            frctr = frctr+1;
            frameTemp = rgb2gray(readFrame(v));
            frameTempRS = imresize(frameTemp,[1080/5 1920/5]);
            vidFrame(:,:,frctr) = frameTempRS(13:end-12,97:end-96);
%             image(vidFrame, 'Parent', currAxes);
%             currAxes.Visible = 'off';
%             pause(1/v.FrameRate);
        end
        
        