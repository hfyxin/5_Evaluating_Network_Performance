%% Generates images contain detected objects, convert to birds eye view. 
% plotting variance oval as well.
%
% remainder of code can be found after 25:00 of video found at:
% https://www.mathworks.com/videos/introduction-to-automated-driving-system-toolbox-1501177087798.html?elqsid=1533134550499&potential_use=Student
clear

% Does data have groundTruth and detection data? y/n = 1/0
groundTruth = 0;     % doesn't work
detected = 0;        % doesn't work
groundTruthWithDetected = 1;   % works
pauseTime = 0;       % pause btw showing images
dispImages = 1;     % display images. 0: only do calculation.
saveImages = 0;     % save images to folder

imgSavePath = '.\results\birds_eye_w_covariance';
resultSavePath = '.\results';
IoU = '0p4';
probErr = 0.95;

% This version only uses rectTPTable.mat
if groundTruth
	% load detections
    % load("C:\Users\benmi\Documents\Thesis\Algorithm Performance\Ground Truth Data\Combined_No_Stationary_Vehicle\Test\testLabels.mat");
    % IoUtruth = testLabelData; % was detections
end

if detected
	% This path is for test data, but not used in later script.
    % load("C:\Users\benmi\Documents\Thesis\Algorithm Performance\YOLOv2\7 Classes Trained on All Data\epoch-step 29250\All Data\IoUdetectedWithConfidence.mat")
	load(".\results\IoUdetectedWithConfidence.mat");
    IoUdetected = IoUdetectedWithConfidence;
end

% get these using IoU_calc_per_frame_v6.m (fixed the FP and FN removal
% problem from v5 by adding FNref and FPref variables
if groundTruthWithDetected
    % These contain ground truth, prediction and ?
    load(".\results\rectFNTable.mat");
    load(".\results\rectFPTable.mat");
    load(".\results\rectTPTable.mat");
end


% location of images
imgPath = '..\All_Combined\Validation';

% Classes to include? y/n = 1/0
CRO = 1;
BIK = 1;
PDU = 1;
PDC = 1;
PDS = 1;
CAR = 0;
BUS = 0;
DGS = 1;
DGL = 1;
SUV = 0;

j=1;
% remove unwanted classes
if CRO == 0
    % IoUtruth.CRO = [];
    IoUdetected.CRO = [];
else
    colour{j} = 'green';
    j=j+1;
end
if BUS == 0
    % IoUtruth.BUS = [];
    IoUdetected.BUS = [];
else    
    colour{j} = 'blue';
    j=j+1;
end
if CAR == 0
    % IoUtruth.CAR = [];
    IoUdetected.CAR = [];
else    
    colour{j} = 'blue';
    j=j+1;
end
if DGS == 0    
    % IoUtruth.DGS = [];
    IoUdetected.DGS = [];
else    
    colour{j} = 'blue';%[0 0.4470 0.7410];
    j=j+1;
end
if DGL == 0      
    % IoUtruth.DGL = [];
    IoUdetected.DGL = [];
else    
    colour{j} = [0.9290 0.6940 0.1250];
    j=j+1;
end
if BIK == 0
    % IoUtruth.BIK = [];
    IoUdetected.BIK = [];
else
    colour{j} = 'black';
    j=j+1;    
end
if PDS == 0
    % IoUtruth.PDS = [];
    IoUdetected.PDS = [];
else    
    colour{j} = 'cyan';
    j=j+1;
end
if PDC == 0
    % IoUtruth.PDC = [];
    IoUdetected.PDC = [];
else
    colour{j} = 'magenta';
    j=j+1;
end
if PDU == 0
    % IoUtruth.PDU = [];
    IoUdetected.PDU = [];
else    
    colour{j} = 'red';
    j=j+1;
end
if SUV == 0
    % IoUtruth.SUV = [];
    IoUdetected.SUV = [];    
end

detections = IoUdetected; % for testing
detGround = rectTPTable;
clear IoUdetected rectTPTable
clear CRO BIK PDU PDC PDS CAR BUS DGS DGL SUV


%% get camera parameters
camParameters = open('.\parameters\cameraParametersJune7Calib.mat');
focalLength = camParameters.cameraParamsJune7Calibration.FocalLength; % [fx, fy] in pixel units
principalPoint = camParameters.cameraParamsJune7Calibration.PrincipalPoint; % [cx, cy] optical center in pixel coordinates
RadialDistortion = camParameters.cameraParamsJune7Calibration.RadialDistortion; %[0 0]
imageSize = [512, 640]; % [nrows, mcols]

camIntrinsics = cameraIntrinsics(focalLength, principalPoint, imageSize, 'RadialDistortion', RadialDistortion);

% Camera tilt angle in degrees
tiltDeg = 8.584;
tilt = tiltDeg * pi/180;

height = 1.67; % mounting height in meters from the ground
pitch  = tiltDeg; % pitch of the camera in degrees
yaw = 0;
roll = 0;
sensorLocation = [-2.235 0.152];

sensor = visionDetectionGenerator('SensorLocation',...
    sensorLocation, 'Height', height,...
    'Pitch', pitch, 'Intrinsics', cameraIntrinsics(...
    focalLength,... %Focal length
    principalPoint,... %principal point
    imageSize,...%image size
    'RadialDistortion',RadialDistortion,...
    'TangentialDistortion',[0 0]),...
    'UpdateInterval',0.1, ... % updated every 0.1s
    'BoundingBoxAccuracy', 5, ...
    'MaxRange', 150); %,...'ActorProfiles', actorProfiles(s)

% used to get the world coordinates from image frame
sensor2 = monoCamera(camIntrinsics,height,'Pitch',pitch, 'sensorLocation', sensorLocation);

% clear them to make workspace cleaner.
clear focalLength principalPoint RadialDistortion imageSize
clear tiltDeg tilt height pitch yaw roll sensorLocation

% since all detections were made on an image of size [700,875] and
% calibration was done for an image of size [512,640], the images need to
% be downsampled back to [512, 640].
scale = 0.73142857;

    
%% Create detection plotter for each class
classes = {detGround.Properties.VariableNames{2:end}};


%% initialize the required plots:
% Use the imageToVehicle() function!!!
% Create extrinsic matrix for predicting object locations 
% rotMatWorld = rotation(pitch, yaw, roll);
% transMatWorld = [sensorLocation(1) sensorLocation(2) height]';
% rotMatCam = rotMatWorld.';
% transMatCam = -rotMatCam*transMatWorld;
% extrinsicMat = [rotMatCam transMatCam; 0 0 0 1];
% camMatrix = cameraMatrix(camParameters.cameraParamsJune7Calibration, rotMatCam, transMatCam)
% camMat2 = [camParameters.cameraParamsJune7Calibration.IntrinsicMatrix zeros(3,1)]*extrinsicMat

%% initialize the required plots:
if dispImages
    % area to display image
    ax1 = axes('Position', [0.02 0 0.55 1]);

    % area to display birds eye plot
    ax2 = axes('Position', [0.62,0.11,0.31,0.82]);
    bep = birdsEyePlot('Parent', ax2,...
        'Xlimits', [0 10],...
        'Ylimits', [-4 4]);
    grid on

    % V shaped area
    covPlot = coverageAreaPlotter(bep,...
        'FaceColor','blue',...
        'EdgeColor','blue');
    % Update coverage area plotter
    plotCoverageArea(covPlot,...
        sensor.SensorLocation,sensor.MaxRange,...
        sensor.Yaw,sensor.FieldOfView(1))

    % create lane marking plotter
    LanePlotter = laneBoundaryPlotter(bep, 'Color', 'red');
    % constant lanes for understanding of object locations
    % driving lane, assuming a lane width of 3m
    lb = parabolicLaneBoundary([-0.00,0.0, 1.5]); 
    rb = parabolicLaneBoundary([-0.0,0.0,-1.5]);
    % adjacent lanes
    lb2 = parabolicLaneBoundary([-0.00,0.0, 4.5]);
    rb2 = parabolicLaneBoundary([-0.0,0.0,-4.5]);
    % Update lanes
    plotLaneBoundary(LanePlotter, [lb lb2, rb rb2])

    % Plot legend
    % legend will not be shown in paper, but for readability in the code.
    for m = 1:length(classes)
        detPlot{m} = detectionPlotter(bep,...
            'MarkerFaceColor',colour{m},...
            'DisplayName', classes{m},...
            'Marker','o');
        detPlotTruth{m} = detectionPlotter(bep,...
            'MarkerFaceColor',colour{m},...
            'DisplayName', classes{m},...
            'Marker','^');
        detPlotVar{m} = detectionPlotter(bep,...
            'MarkerFaceColor','non',...
            'DisplayName', classes{m},...
            'Marker','.');
        legend
    end

    truthPlot = outlinePlotter(bep);
    % set figure size
    set(gcf, 'Position', [10 10 1200 600])
end

%% parse through all images and plot their detections with bounding boxes
% and world coordinates on the plot.

fprintf("Processing detection table with %d rows:\n", size(detGround,1));
for k = 1:size(detGround, 1)  % each image
% size(detGround, 1)-4% 370:390%
%     img = imread([imgPath,'\',detections.imageFilename{k}])*scale;
%     imgResize = img*scale;
%     imshow(img,'Parent',ax1)

    % init
    boundingBox = cell(1,length(classes));
    boundingBoxTruth = cell(1,length(classes));
    % get bbox for current image
    for m = 1:length(classes)     % each class (column)

        % get the bounding Box of the prediction
        for h = 1:size(detGround{k,:}{m+1},1)   % each detection of same class
            if size(detGround{k,:}{m+1}{h,3},1)>0  % detection found
                % bounding box dimension
                boundingBox{m} = [boundingBox{m} ; detGround{k,:}{m+1}{h,3}{1}(:,1:4)*scale];%detections{k,:}{m+1}(:,1:4)
                % minimum bounding box value can be is 1
                boundingBox{m}(boundingBox{m}<1)=1;
            end
        end
        
        % get the bounding Box of the ground truth
        for g = 1:size(detGround{k,:}{m+1},1)
            if size(detGround{k,:}{m+1}{g,1},1)>0
                count = 1;
                for t = 1:size(detGround{k,:}{m+1}{g,1}{1},1)
                    if max(detGround{k,:}{m+1}{g,1}{1}(t,1:4))==0

                        boundingBoxTruth{m} = [boundingBoxTruth{m} ; detGround{k,:}{m+1}{g,1}{1}(t-count,1:4)*scale];%detections{k,:}{m+1}(:,1:4)
                        % minimum bounding box value can be is 1
                        boundingBoxTruth{m}(boundingBoxTruth{m}<1)=1;
                        count = count+1;
                    else
                        boundingBoxTruth{m} = [boundingBoxTruth{m} ; detGround{k,:}{m+1}{g,1}{1}(t,1:4)*scale];
                        boundingBoxTruth{m}(boundingBoxTruth{m}<1)=1;
                    end
                end
            end
        end
        
    end
    
    % find the location of the object, i.e. the bottom center point.
    objLocation = cell(1,length(classes));
    objLocationTruth = cell(1,length(classes));    
    for m = 1:length(classes)
        for n = 1:size(boundingBox{m},1)
            objLocation{m}(n, :) = [boundingBox{m}(n, 1)+boundingBox{m}(n, 3)/2 , boundingBox{m}(n, 2)+boundingBox{m}(n, 4)];
        end
        for p = 1:size(boundingBoxTruth{m},1)
            objLocationTruth{m}(p, :) = [boundingBoxTruth{m}(p, 1)+boundingBoxTruth{m}(p, 3)/2 , boundingBoxTruth{m}(p, 2)+boundingBoxTruth{m}(p, 4)];
        end
    end
    
    % handle out-of-image location points
    for m = 1:length(classes)
        if length(objLocation{m})>0
            objLocationCol{m} = objLocation{m}(:,2);
            objLocationCol{m}(objLocationCol{m}>512)=512;
            objLocation{m}(:,2) = objLocationCol{m};
            objLocationRow{m} = objLocation{m}(:,2);
            objLocationRow{m}(objLocationRow{m}>640)=640;
            objLocation{m}(:,2) = objLocationRow{m};
            
            objLocationColTruth{m} = objLocationTruth{m}(:,2);
            objLocationColTruth{m}(objLocationColTruth{m}>512)=512;
            objLocationTruth{m}(:,2) = objLocationColTruth{m};
            objLocationRowTruth{m} = objLocationTruth{m}(:,2);
            objLocationRowTruth{m}(objLocationRowTruth{m}>640)=640;
            objLocationTruth{m}(:,2) = objLocationRowTruth{m};
        end
    end
    
    % Map from pixel to world coordinate
    % Add world coordinate location of ground truth and detected object for
    % each true positive (didn't want to figure out indexing so making
    % detected world-coordinate and ground truth world-coordinate pairs instead
    for m = 1:length(classes)
        if size(objLocation{m}, 1)>0
            detectionsWorld{m} = imageToVehicle(sensor2, objLocation{m});
            detectionsWorldTruth{m} = imageToVehicle(sensor2, objLocationTruth{m});
            % if any detections come back as > 500m, make their values 500m
            % if they come back negative, make their values 500m 
            % test(test>999)=999
            detComb{m} = {detectionsWorldTruth{m}, detectionsWorld{m}};
            detPixComb{m} = {objLocationTruth{m}, objLocation{m}};
            
            % For debugging negative detection problem.
			% Exception found here: row-285/frame_1899, possibly caused by
			% a pedestrian on topleft being too far away (uneven ground). -Elliot
            if detectionsWorld{m}(:,1)<0
                fprintf("<Negative detection row-%d>, ", k);
            end
        else
            detectionsWorld{m} = [];
            detectionsWorldTruth{m} = [];
            detComb{m} = {[]};
            detPixComb{m} = {[]};
        end
    end
    
    % IMPORTANT. Store location coordinate pairs into a table.
	% the pair looks like: [ground truth, detected]
    if k==1 || exist('WCPairTable')==0
        WCPairTable = array2table(detComb,'VariableNames',classes);
        pixPairTable = array2table(detPixComb,'VariableNames',classes);
    else
        WCPairTableAdd = array2table(detComb,'VariableNames',classes);
        pixPairTableAdd = array2table(detPixComb,'VariableNames',classes);
        
        WCPairTable = [WCPairTable; WCPairTableAdd];
        pixPairTable = [pixPairTable; pixPairTableAdd];
    end

    % show image and plot, configurable
    if dispImages
        % Plot the video frame and bounding boxes for the objects as well as the
        % BEP
        img = imread([imgPath,'\',detGround.imageFilename{k}]);% detections.imageFilename{k}]);
        frameAnnotated = imresize(img,scale);

        i=0;
        % clear all detections shown on detPlot from previous frame
        for m=1:length(classes)  
            clearData(detPlot{m})
            clearData(detPlotTruth{m})
            clearData(detPlotVar{m})
        end

        for m = 1:length(classes)    
            if size(objLocation{m},1) > 0
                % draw bounding boxes and location point on image
                frameAnnotated = insertShape(frameAnnotated, 'Rectangle', ...
                    boundingBoxTruth{m}, 'Color', 'yellow', 'LineWidth', 3);
                frameAnnotated = insertShape(frameAnnotated, 'Rectangle', ...
                    boundingBox{m}, 'Color', colour{m}, 'LineWidth', 3);
                    %[leftmost x,top y,width, height]
                % textCoord = [num2str(objLocation{m}(1),3), ', ', ...
                %             num2str(objLocation{m}(2),3)];
                % frameAnnotated = insertText(frameAnnotated, objLocation{m}+[10,10], ...
                %    textCoord, 'BoxColor','white', 'BoxOpacity', 0.9);
                    %colour{m});%[leftmost x,top y,width, height]
                    
                % draw confidence oval, point by point
                for t = 1:size(objLocation{m},1)
                    [corr_loc, radius] = convertRawPCSToWCSFcn(objLocation{m}(t,:), probErr);
                    x0 = corr_loc(1); y0 = corr_loc(2);
                    t=-pi:0.2:pi;
                    x=x0+radius(1)*cos(t);
                    y=y0+radius(2)*sin(t);
                    % draw oval
                    % frameAnnotated = insertShape(frameAnnotated, 'FilledCircle', ...
                    % [x',y',1*ones(length(x'),1)],'Color','green','Opacity',1);
                    
                    % draw oval
                    frameAnnotated = insertMarker(frameAnnotated, [x',y'], ...
                        '+', 'Color', 'green', 'size', 1);
                    % draw center
                    frameAnnotated = insertShape(frameAnnotated, 'FilledCircle', ...
                        [x0,y0,3],'Color','green','Opacity',1, 'LineWidth', 1);
                    % draw text
                    textCoord = [num2str(corr_loc(1),3),', ',num2str(corr_loc(2),3)];
                    frameAnnotated = insertText(frameAnnotated, corr_loc +[10,10], ...
                        textCoord, 'BoxColor','white', 'BoxOpacity', 0.8);
                    %colour{m});%[leftmost x,top y,width, height]
                    
                    
                end
            else
            end
            % show image
            im = imshow(frameAnnotated, 'Parent', ax1);
            
            % plots detection on bird's eye plot
            if length(detectionsWorld{m})>0
                plotDetection(detPlot{m}, detectionsWorld{m})
                plotDetection(detPlotTruth{m}, detectionsWorldTruth{m})
                % draw confidence oval (on world coordinate)
                for t = 1:size(objLocation{m},1)     
                    [corr_loc, radius] = convertRawPCSToWCSFcn(objLocation{m}(t,:), probErr);
                    x0 = corr_loc(1); y0 = corr_loc(2);
                    t=-pi:0.2:pi;
                    x=x0+radius(1)*cos(t);
                    y=y0+radius(2)*sin(t);
                    % remove out-of-image points
                    y = y(x<640);
                    x = x(x<640);
                    x = x(y<512);
                    y = y(y<512);
                    y = y(x>1);
                    x = x(x>1);
                    x = x(y>1);
                    y = y(y>1);
                    % pixel to world conversion
                    varianceWorld = imageToVehicle(sensor2, [x',y']);
                    plotDetection(detPlotVar{m}, varianceWorld)
                end
                hold on
                grid on
            end
        end

        if saveImages == 1
            if ~exist(imgSavePath, 'dir')
                mkdir(imgSavePath)
            end
            saveas(gcf, [imgSavePath, '\',detGround.imageFilename{k}, '_IoU_', IoU,'.jpg'])
        end

        % Also possible to plot tracks using trackPlotter(bep). See:
        % https://www.mathworks.com/help/driving/examples/visualize-sensor-coverage-detections-and-tracks.html
        pause(pauseTime)    
    end
    
	% display progress
	if rem(k, 100) == 1
		fprintf("%d, ", k);
	end
end
fprintf("done.\n");

% May consider save pixPairTable and WCPairTable
save([resultSavePath, '\', 'pixPairTable'], 'pixPairTable')
save([resultSavePath, '\', 'WCPairTable'], 'WCPairTable')
fprintf("Coordinate values saved to: \n\tpixPairTable.mat \n\tWCPairTable.mat\n\n")
