%% Plot distribution curves for a grid of XY range (in pixel).
%
% Created by Elliott, based on Ben's files:
%   - sortPixelErrorsXdistributionWithX_full.m
%   - sortPixelErrorsYdistributionWithX_full.m
% 
% The original image of 640 * 512 is divided into 2 * 3 grid (range):
%               +-----------------------------+   --> (x direction)
%    the image: | range 1 | range 2 | range 3 |
%               +---------+---------+---------+
%               | range 4 | range 5 | range 6 |
% (y direction) +-----------------------------+ 
%
% This script is to show that horizontal grids have similar error distribution.
%

%% Load coordinate pairs
% pair_value looks like: [ground truth, detected]
clear
load('.\results\WCPairTable.mat')     % WCPairTable, WCS: world coordinate system
load('.\results\pixPairTable.mat')    % pixPairTable, 

%% Pre-calculated values to divide the image
range_y = [1,324,512]';
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
% pre-set text for all ranges (for plotting)
range_txt = cell(n_ranges, 1);
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
			x_pix = pred(i,1); % x coordinate, use detected value
			y_pix = pred(i,2); % y coordinate
			% find range (grid) index
			idx_x = find(range_x > x_pix, 1) - 1;
			idx_y = find(range_y > y_pix, 1) - 1;
			idx = (idx_y - 1) * grid_c + idx_x;  % e.g.: id(1,2) is converted to 2
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



