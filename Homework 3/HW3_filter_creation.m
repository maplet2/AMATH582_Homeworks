clear; close all; clc;

% Load Data
file_name = 'cam1_1.mat';
camm = load(file_name);
    
% Convert the video data into a matrix of numeric values
% Note that cell2mat will return a 4D uint matrix where the
%   there is a 480 x 640 pixel video that has 226 (in the 
%   cam1_1.mat case) frames and then 3 other matrices that
%   represents the RGB value of the picture
cam = cell2mat(struct2cell(camm));
cam_dim = size(cam);
imsize = cam_dim(1:2);
imx = imsize(2);
imy = imsize(1);
time = cam_dim(4);

% % Play first k frames
% k = 50;
% for n = 1:k
%     frame = cam(:,:,:,n);
%     imshow(uint8(frame))
%     hold on 
%     pause(0.1)
% end

%%
% Initialize Filter Bounds

% Camera 1 Filter Bounds
fxmin = 250;
fxmax = 435;
fymin = 200;
fymax = 400;

% % Camera 2-1, 2-2 Filter Bounds
% fxmin = 200;
% fxmax = 400;
% fymin = 50;
% fymax = 400;

% % Camera 2-3 Filter Bounds
% fxmin = 200;
% fxmax = 400;
% fymin = 150;
% fymax = 400;

% % Camera 2-4 Filter Bounds
% fxmin = 200;
% fxmax = 390;
% fymin = 95;
% fymax = 350;

% % Camera 3 Filter Bounds
% fxmin = 250;
% fxmax = imy;
% fymin = 150;
% fymax = 350;

% Create filter
filter = zeros(imy, imx); 
filter(fymin:fymax, fxmin:fxmax) = ones((fymax-fymax)+1, (fymax-fymax)+1);
minyy = zeros(1,time);      % Store object's y location
minxx = zeros(1,time);      % Store object's x location
threshold = 30;

% Note the target choice for each different camera is as follows:
%   Camera 1: The Reddish-Hue of the Flashlight in all cases
%   Camera 2: 2-1 White, 2-2 Blue-White, 2-3 and 2-4 Reddish Hue
%   Camera 3: 3-1 3-2 3-3 White
for n = 1:time
    % Change frame values into doubles
    frame = double((cam(:,:,:,n)));
    
    % Show Filter on Video as necessary
    frame(:,:,1) = frame(:,:,1).*filter;
    frame(:,:,2) = frame(:,:,2).*filter;
    frame(:,:,3) = frame(:,:,3).*filter;
    imshow(uint8(frame));
    pause(0.1)

    % Choose Target
    target = [223 154 151];       % Reddish Hue of Flashlight
    % target = [255 255 255];        % White
    % target = [243 254 254];        % Blue-White emitted from flashlight
    
    % Find Location in Filtered Frame that's the closest
    %   to our Target color
    dist = colordiff(frame,filter, target); 
%     dist2 = colordiff(r,g,b,[255 255 255]);
%     if min(dist(:)) > threshold
%         dist = dist2;
%     end
    [min_val, min_ind] = min(dist(:));
    [miny, minx] = ind2sub(imsize, min_ind);
    
    % Record x and y location of 'closest pixel'
    minyy(n) = miny;
    minxx(n) = minx;           
end

% Save x and y locations into a .mat file
% save('camera1_1.mat', 'minxx','minyy')
%%

% The following code will check path of recorded 
%   against the original video 
for n = 1:time
    frame = cam(:,:,:,n);
    framebw = rgb2gray(frame);
    imshow(uint8(frame))
    hold on 
    plot(minxx(n), minyy(n), 'r.', 'markersize',20)
    pause(0.1)
end

% Plots the path of the oscillating object.
figure(2) 
plot(minxx, imy-minyy, 'r.', 'markersize', 20)
xlim([0 imx])
ylim([0 imy])


% Receives a frame and target color and then calculates
%   the difference between the target color and each pixel 
%   of the frame and returns a matrix of differences
function dist = colordiff(frame, filter, target)
    r = frame(:,:,1).*filter;
    g = frame(:,:,2).*filter;
    b = frame(:,:,3).*filter;
    dist = sqrt((r-target(1)).^2 + (g-target(2)).^2 + (b-target(3)).^2);
end