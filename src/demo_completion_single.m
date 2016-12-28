function [] = demo_completion_single(frags, imgNum, params)
% get two inducers from the user and complete a single curve between them


close all

imgName = params.imgNames{imgNum};
baseName = imgName(1:end-4);
img = imread([params.imgsFolder imgName]);
imgRGB = img;

imshow(imgRGB);

end

