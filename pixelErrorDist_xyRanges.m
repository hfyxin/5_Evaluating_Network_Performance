%% Plot distribution curves for a grid of XY range (in pixel).
%
% Suggest checking pixelErrorDist_yRanges.m before this one. It's easier.
% Created by Elliott, based on Ben's files:
%   - sortPixelErrorsXdistributionWithX_full.m
%   - sortPixelErrorsYdistributionWithX_full.m
% The change is verified by comparing the result from Ben's code.

% The original image of 640 * 512 is divided into 2 * 3 grid (range):
%               +--------------------------+   --> (x direction)
%    the image: | range1 | range2 | range3 |
%               +--------+--------+--------+
%               | range4 | range5 | range6 |
% (y direction) +--------------------------+ 
%
% This script is to show if horizontal grids have similar error distribution.
%

%% Load coordinate pairs
% pair_value looks like: [ground truth, detected]
clear
load('.\results\WCPairTable.mat')     % WCPairTable, WCS: world coordinate system
load('.\results\pixPairTable.mat')    % pixPairTable, 

%% Pre-calculated values to divide the image
% range_y = [1,324,512]';
range_y = [1,253,512]';
range_x = [1,214,427,640]';

%% init some parameters
range_y(end) = range_y(end) + 1;   % to incorperate the largest pixel value
range_x(end) = range_x(end) + 1;   % to incorperate the largest pixel value
grid_r = length(range_y) - 1;
grid_c = length(range_x) - 1;
hist_grid = [grid_r,grid_c];       % the [row, col] of multi plots (histogram)
n_ranges = grid_r * grid_c;
pixErrY = cell(n_ranges, 1);      % y pixel error array
pixErrX = cell(n_ranges, 1);      % x pixel error array
pixErrDistY = cell(n_ranges, 1);  % distribution function for y
pixErrDistX = cell(n_ranges, 1);  % distribution function for x

% pre-set label text for all ranges (for plotting)
range_txt = cell(n_ranges, 1);   % cell, contains x, y range in str format.
for r = 1:grid_r
    for c = 1:grid_c
        i = (r-1)*grid_c + c;
        y_txt = ['[', num2str(range_y(r)), ',', num2str(range_y(r+1)), ')'];
        x_txt = ['[', num2str(range_x(c)), ',', num2str(range_x(c+1)), ')'];
        range_txt{i} = {x_txt; y_txt};
    end
end


%% record x and y pixel error values within each range
for row = 1:height(pixPairTable)
    for col = 1:width(pixPairTable)
        % get the table cell
        detection = pixPairTable{row, col}{1};  % cell array, don't forget the last {1}
        [cell_r, cell_c] = size(detection);

        % Possible values of detection cell array:
        %                     | cell size | detection[1] | detection[2] |
        %                     |   r, c    | Ground Truth |  Predicted   |
        % --------------------+-----------+--------------+--------------+
        %  Nothing detected   |   1 x 1   |     []       |      null    |  
        % One object detected |   1 x 2   |  1 x 2 array |  1 x 2 array |
        % Multiple detections |   1 x 2   |  n x 2 array |  n x 2 array |
        
        % skip this cell if contains nothing
        if cell_c == 1
            continue
        end
        
        % detection found, extract ground truth and detection
        gt = detection{1};     % n x 2 array
        pred = detection{2};   % n x 2 array
        for i = 1:size(gt,1)   % iterate through i rows
            y_pix = pred(i,2); % y coordinate, use detected value
            x_pix = pred(i,1); % x coordinate, use detected value
            % find grid index, then range index
            idx_x = find(range_x > x_pix, 1) - 1;
            idx_y = find(range_y > y_pix, 1) - 1;
            idx = (idx_y - 1) * grid_c + idx_x;  % e.g.: id(1,2) is converted to 2
            % record pixel error
            x_err = pred(i,1) - gt(i,1);         % pixel error on x direction
            y_err = pred(i,2) - gt(i,2);         % pixel error on y direction
            pixErrX{idx} = [pixErrX{idx}; x_err];
            pixErrY{idx} = [pixErrY{idx}; y_err];
        end
    end
end


%% Plot Y direction pixel error histograms, separate plot
figure
n_row = hist_grid(1); n_col = hist_grid(2);  % number of plots on each row, col
sgtitle('Pixel error distribution on Y direction') % subplot grid title
for i = 1:n_ranges
    subplot(n_row,n_col,i);
    histogram(pixErrY{i});
    xlabel('pixel error');
    ylabel('# of detections');
    % title example: range[1,161)
    title_1 = ['Range: X=', range_txt{i}{1}];
    title_2 = ['Y=', range_txt{i}{2}];
    title({title_1; title_2});
    ax = gca;
    ax.YGrid = 'on';
end


%% Plot X direction pixel error histograms, separate plot
figure
n_row = hist_grid(1); n_col = hist_grid(2);  % number of plots on each row, col
sgtitle('Pixel error distribution on X direction') % subplot grid title
for i = 1:n_ranges
    subplot(n_row,n_col,i);
    histogram(pixErrX{i});
    xlabel('pixel error');
    ylabel('# of detections');
    % title example: range[1,161)
    title_1 = ['Range: X=', range_txt{i}{1}];
    title_2 = ['Y=', range_txt{i}{2}];
    title({title_1; title_2});
    ax = gca;
    ax.YGrid = 'on';
end


%% Fit pixel error values in each range to probability function (PDF)
for i = 1:n_ranges
    pixErrDistY{i} = fitdist(pixErrY{i}, 'Normal');
    pixErrDistX{i} = fitdist(pixErrX{i}, 'Normal');
end


%% find z score. Need to compare both horizontally and vertically
% prepare some variables
n_objects = zeros(n_ranges, 1);    % number of objects in each range
mu_y = zeros(n_ranges, 1);
mu_x = zeros(n_ranges, 1);
sigma_y = zeros(n_ranges, 1);
sigma_x = zeros(n_ranges, 1);
for i = 1:n_ranges
    n_objects(i) = length(pixErrY{i});
    mu_y(i) = pixErrDistY{i}.mu;
    mu_x(i) = pixErrDistX{i}.mu;
    sigma_y(i) = pixErrDistY{i}.sigma;
    sigma_x(i) = pixErrDistX{i}.sigma;
end

% z score init
z_score_hori_y = zeros(n_ranges, 1);     % horizontally, w.r.t the next grid
z_score_hori_x = zeros(n_ranges, 1);     % horizontally, w.r.t the next grid
z_score_vert_y = zeros(n_ranges, 1);     % vertically, w.r.t the next grid
z_score_vert_x = zeros(n_ranges, 1);     % vertically, w.r.t the next grid

for i = 1:n_ranges
    % skip the last col, when calculating horizontal z score
    if rem(i, grid_c) == 0
        z_score_hori_y(i) = nan;
        z_score_hori_x(i) = nan;
    else
        % calculate horizontally adjacent grid (e.g., range1 vs range2)
        z_score_hori_y(i) = (mu_y(i) - mu_y(i+1)) / ...
            sqrt(sigma_y(i)^2/n_objects(i) + sigma_y(i+1)^2/n_objects(i+1)); 
        z_score_hori_x(i) = (mu_x(i) - mu_x(i+1)) / ...
            sqrt(sigma_x(i)^2/n_objects(i) + sigma_x(i+1)^2/n_objects(i+1)); 
    end
    % skip the last row, when calculating vertical z score
    if i > n_ranges - grid_c
        z_score_vert_y(i) = nan;
        z_score_vert_x(i) = nan;
    else
        % calculate vertically adjacent grid (e.g., range1 vs range4)
        z_score_vert_y(i) = (mu_y(i) - mu_y(i+grid_c)) / ...
            sqrt(sigma_y(i)^2/n_objects(i) + sigma_y(i+grid_c)^2/n_objects(i+grid_c)); 
        z_score_vert_x(i) = (mu_x(i) - mu_x(i+grid_c)) / ...
            sqrt(sigma_x(i)^2/n_objects(i) + sigma_x(i+grid_c)^2/n_objects(i+grid_c)); 
    end
end


%% Summarize what we've got so far into a table
% note the var names begins with Capital to look better.
Range = (1:n_ranges)';    % range ID
XRange = cell(n_ranges,1);  % str text cell
YRange = cell(n_ranges,1);  % str text cell
for i = 1:n_ranges
    XRange{i} = range_txt{i}{1};
    YRange{i} = range_txt{i}{2};
end
% the summary table
t_dist = table(Range,XRange,YRange,n_objects,...
                mu_y,sigma_y,z_score_hori_y,z_score_vert_y,...
                mu_x,sigma_x,z_score_hori_x,z_score_vert_x);
t_dist


%% Plot Y direction pixel error distribution curve, single plot
figure
x_ = -40:0.1:60;
legend_val = cell(1,n_ranges);
for i = 1:n_ranges
    % plot PDF
    y_ = pdf(pixErrDistY{i}, x_);
    plot(x_, y_, 'LineWidth', 2);
    legend_txt = ['Range: X=', range_txt{i}{1}, ' Y=', range_txt{i}{2}];
    legend_val{i} = legend_txt;
    hold on
    % plot mu
    ver_line = pixErrDistY{i}.mu;
    % xline(ver_line, '--k');   % only in MATLAB 2018b or higher
end
title('Fitted PDF for errors on Y direction');
legend(legend_val,'Location', 'NorthEast');
xlabel('pixel error');
ylabel('probability');
grid on

%% Plot X direction pixel error distribution curve, single plot
figure
x_ = -40:0.1:60;
legend_val = cell(1,n_ranges);
for i = 1:n_ranges
    % plot PDF
    y_ = pdf(pixErrDistX{i}, x_);
    plot(x_, y_, 'LineWidth', 2);
    legend_txt = ['Range: X=', range_txt{i}{1}, ' Y=', range_txt{i}{2}];
    legend_val{i} = legend_txt;
    hold on
    % plot mu
    ver_line = pixErrDistY{i}.mu;
    % xline(ver_line, '--k');   % only in MATLAB 2018b or higher
end
title('Fitted PDF for errors on X direction');
legend(legend_val,'Location', 'NorthEast');
xlabel('pixel error');
grid on

%% remove some temp vars
clearvars -except pixPairTable WCPairTable range_y range_x hist_grid n_ranges ...
                  pixErrX pixErrY pixErrDistX pixErrDistY t_dist range_txt
