%% Plot distribution curves for each range in Y direction (in pixel).
%
% Created by Elliot Huangfu, based on Ben's files, but with some addition.
%   - sortPixelErrorsXdistYpixels_full.m
%   - sortPixelErrorsYdistribution_full.m
% The change is verified by comparing the result from Ben's code.
% However, some table does not perfect match the paper draft (but close).
%
% This script is to show that error distributions vary by distance (y range).
%

%% Load coordinate pairs
% pair_value looks like: [ground truth, detected]
clear
load('.\results\WCPairTable.mat')     % WCPairTable, WCS: world coordinate system
load('.\results\pixPairTable.mat')    % pixPairTable, 

%% Pre-calculated values to divide the image, don't forget it starts from 1
% original pixel range, calculated by WCS range [1,5,10,15, ...,55,60,75,90,inf]
range_y_original = [1 161 165 168 172 175 178 182 188 195 206 223 253 324 512]';   % 14 ranges
% the range from Ben's code sortPixelErrorsYdistribution_full.m
range_y_pick1 = [1 161 195 324 512]';      % 4 ranges
% the range Elliot observed
range_y_pick2 = [1 178 195 223 324 512]';  % 5 ranges
% the range random picked
range_y_pick3 = [1 195 512]';

% SELECT YOUR CHOICE HERE
ranges = range_y_pick2;      % which range setting to process
hist_grid = [2,3];           % the [row, col] of multi plots (histogram)

%% init some parameters
ranges(end) = ranges(end) + 1;   % to incorperate the largest pixel value
n_ranges = length(ranges) - 1;   % number of ranges
pixErrY = cell(n_ranges, 1);     % y pixel error array
pixErrX = cell(n_ranges, 1);     % x pixel error array
pixErrDistY = cell(n_ranges, 1);  % distribution function for y
pixErrDistX = cell(n_ranges, 1);  % distribution function for x

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
            % find range index, neat!
            idx = find(ranges > y_pix, 1) - 1;
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
    title(['range [', num2str(ranges(i)), ',', num2str(ranges(i+1)), ')']);
    % axis([-40 40 0 75])
    ax = gca;
    ax.YGrid = 'on';
end


%% Fit pixel error values in each range to probability function (PDF)
for i = 1:n_ranges
    pixErrDistY{i} = fitdist(pixErrY{i}, 'Normal');
    pixErrDistX{i} = fitdist(pixErrX{i}, 'Normal');
end


%% find z score
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
z_score_y = zeros(n_ranges, 1);     % w.r.t the next range
z_score_x = zeros(n_ranges, 1);     % w.r.t the next range
for i = 1:n_ranges
    % skip the last calculation
    if i == n_ranges
        z_score_y(i) = nan;
        z_score_x(i) = nan;
        continue
    end
    % calculation
    z_score_y(i) = (mu_y(i) - mu_y(i+1)) / ...
        sqrt(sigma_y(i)^2/n_objects(i) + sigma_y(i+1)^2/n_objects(i+1)); 
    z_score_x(i) = (mu_x(i) - mu_x(i+1)) / ...
        sqrt(sigma_x(i)^2/n_objects(i) + sigma_x(i+1)^2/n_objects(i+1)); 
end
    

%% Summarize what we've got so far into a table
% note the var names begins with Capital to look better.
Range = (1:n_ranges)';    % range ID
PixelRangeStart = ranges(1:end-1);  % inclusive
PixelRangeEnd = ranges(2:end);    % exclusive
% the table used in paper
t_dist = table(Range,PixelRangeStart,PixelRangeEnd,n_objects,...
                mu_y,sigma_y,z_score_y, mu_x,sigma_x,z_score_x);
t_dist


%% Plot Y direction pixel error distribution curve, single plot
figure
x_ = -40:0.1:60;
legend_val = cell(1,n_ranges);
for i = 1:n_ranges
    % plot PDF
    y_ = pdf(pixErrDistY{i}, x_);
    plot(x_, y_, 'LineWidth', 2);
    legend_val{i} = ['range [', num2str(ranges(i)), ',', num2str(ranges(i+1)), ')'];
    hold on
    % plot mu
    ver_line = pixErrDistY{i}.mu;
    % xline(ver_line, '--k');   % only in MATLAB 2018b or higher
end
title('Fitted PDF for errors on Y direction');
legend(legend_val,'Location', 'NorthEast');
xlabel('pixel error');
grid on

%% remove some temp vars
clearvars -except pixPairTable WCPairTable ranges hist_grid n_ranges ...
                  pixErrX pixErrY pixErrDistX pixErrDistY t_dist