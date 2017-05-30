function imOut = imgLandoltC(varargin)
% Build a Landolt C image
% 
% 
%%

p = inputParser;
p.addParameter('orientation','up',@ischar);
p.addParameter('outersize',100,@isnumeric);
p.addParameter('innersize',50,@isnumeric);
p.addParameter('gapsize',10,@isnumeric);
p.parse(varargin{:});

orientation = p.Results.orientation;
outerSize = p.Results.outersize;
innerSize = p.Results.innersize;
gapSize = p.Results.gapsize;

th = (0:.01:1)'*2*pi;
xh = outerSize/2*cos(th);
yh = outerSize/2*sin(th);

xc = innerSize/2*cos(th);
yc = innerSize/2*sin(th);

imL = zeros(outerSize+1,outerSize+1);
[mr,mc] = meshgrid([-outerSize/2:outerSize/2]);
radv = sqrt(mr.^2 + mc.^2);
imL(radv>innerSize/2 & radv<outerSize/2) = 1;

cgap = gapSize;

switch orientation
    case 'left'
        imL(outerSize/2-gapSize/2:outerSize/2+gapSize/2,1:outerSize/2) = 0;
    case 'right'
        imL(outerSize/2-gapSize/2:outerSize/2+gapSize/2,outerSize/2:end) = 0;
    case 'up'
        imL(1:outerSize/2,outerSize/2-gapSize/2:outerSize/2+gapSize/2) = 0;
    case 'down'
        imL(outerSize/2:end,outerSize/2-gapSize/2:outerSize/2+gapSize/2) = 0;
end


imOut = imL; 

% cNoise = 10;
% imOutBig = zeros(101+cNoise,101+cNoise);
% xs = ceil(cNoise*(-.5+rand(1,1)));
% ys = ceil(cNoise*(-.5+rand(1,1)));
% imOutBig(xs+cNoise/2+[1:101],ys+cNoise/2+[1:101]) = imL;
% imOut(1:101,1:101) = imOutBig(cNoise/2+[1:101],cNoise/2+[1:101]);
% figure; imagesc(imOut);