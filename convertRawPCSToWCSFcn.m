function [pix_corr, radius] = convertRawPCSToWCSFcn(pixel, probErr)
% Convert detected pixel coords to estimated pixel coords w/ oval radius.
% Created by Ben. Modified by Elliot Huangfu to make it configurable.

% original range by Ben
range1 = 195;
range2 = 206;
range3 = 223;
range4 = 253;
range5 = 324;
range6 = 512;

% a better way to do it.
% ranges = [1,195,206,223,253,324,512];   % origianl range by Ben
ranges = [1,178,195,223,324,512];       % revised range by Elliot
ranges(end) = ranges(end) + 1;

% mean, variance of the error, for x and y error respectively.
distributions = struct;
distributions.mu = [-0.57883, -1.3636
                    -0.36435, -2.1413
                    -0.92757, -0.99214
                    -1.0114,  1.536
                    -1.8435,  11.756];
distributions.sigma = [3.7648, 4.0471
                       3.5759, 4.8228
                       3.5794, 5.9277
                       4.8948, 9.7961
                       6.3796, 19.122];

% Out of range error
if pixel(2) >= ranges(end) || pixel(2) < ranges(1)
    error(['pixel y coord out of range: ', num2str(pixel(2))]);
end

% find the range that y coord falls into, Elliot
idx = find(ranges > pixel(2), 1) - 1;
err = distributions.mu(idx,:);
sigma = distributions.sigma(idx,:);

% find the range that y coord falls into, Ben
% if pixel(2)<range1
%     err = [-.345569,-1.9];
%     sigma = [3.6243, 4.60386];
% elseif pixel(2)>=range1 && pixel(2) < range2
%     err = [-.345569,-0.820233];
%     sigma = [3.6243, 5.67147];
% elseif pixel(2)>=range2 && pixel(2) < range3
%     err = [-1.31711,0.316178];
%     sigma = [4.05014, 8.89495];
% elseif pixel(2)>=range3 && pixel(2) < range4
%     err = [-.187478,0.316178];
%     sigma = [4.51462, 8.89495];
% elseif pixel(2)>=range4 && pixel(2) < range5
%     err = [-.187478,0.316178];
%     sigma = [5.90608, 8.89495];
% elseif pixel(2)>=range5 && pixel(2) < range6
%     err = [-1.769,10.2578];
%     sigma = [5.90608, 17.5366];
% end

% choose the probable error, confidence level
if probErr == 0.5
    std_dev = 0.6745;
elseif probErr == 0.95
    std_dev = 1.96;
end

radius = sigma*std_dev;
pix_corr = pixel-err;

% px_min = pix_corr(1)-std_dev*sigma(1);
% px_max = pix_corr(1)+std_dev*sigma(1);
% py_min = pix_corr(2)-std_dev*sigma(2);
% py_max = pix_corr(2)+std_dev*sigma(2);
% 
%     
% camParameters = open('C:\Users\benmi\Documents\Thesis\Matlab\Calibration Parameters\cameraParametersJune7Calib.mat');
% focalLength = camParameters.cameraParamsJune7Calibration.FocalLength; % [fx, fy] in pixel units
% principalPoint = camParameters.cameraParamsJune7Calibration.PrincipalPoint; % [cx, cy] optical center in pixel coordinates
% RadialDistortion = camParameters.cameraParamsJune7Calibration.RadialDistortion; %[0 0]
% imageSize = [512, 640]; % [nrows, mcols]
% camIntrinsics = cameraIntrinsics(focalLength, principalPoint, imageSize, 'RadialDistortion', RadialDistortion);
% 
% % Camera tilt angle in degrees
% tiltDeg = 8.584;
% tilt = tiltDeg * pi/180;
% 
% height = 1.67; % mounting height in meters from the ground
% pitch  = tiltDeg; % pitch of the camera in degrees
% yaw = 0;
% roll = 0;
% sensorLocation = [-2.235 0.152];
% sensor2 = monoCamera(camIntrinsics,height,'Pitch',pitch, 'sensorLocation', sensorLocation);
% % since all detections were made on an image of size [700,875] and
% % calibration was done for an image of size [512,640], the images need to
% % be downsampled back to [512, 640].
% scale = 0.73142857;
% 
% xy_base = imageToVehicle(sensor2, pix_corr)
% x_max = imageToVehicle(sensor2, [pix_corr(1), py_min])
% x_min = imageToVehicle(sensor2, [pix_corr(1), py_max])
% y_min = imageToVehicle(sensor2, [px_max, pix_corr(2)])
% y_max = imageToVehicle(sensor2, [px_min, pix_corr(2)])
end
