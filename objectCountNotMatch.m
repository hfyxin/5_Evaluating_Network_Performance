%% Proof that Ben's code has bug processing the true positives.
% script DisplayBirdsEyeGTruthAndDet_v7.m
%   - takes in True Positive results: rectTPTable.mat
%   - outputs the location (bottom center) of each bounding box
% 
% This script shows that the number of objects don't match.
% The exact reason has not been found, but for now will rely on
% rectTPTable.m
%

%% The two tables
clear
load(".\results\rectTPTable.mat");   % only care about True Positives
rectTPTableClean = rectTPTable(:,2:end);  % drop the first column
load('.\results\pixPairTable.mat');  % pixPairTable

%% check the labels of 2 tables match
labelTPTable = rectTPTableClean.Properties.VariableNames;      % col names
labelPairTable = pixPairTable.Properties.VariableNames;   % col names
assert(strcmp([labelTPTable{:}],[labelPairTable{:}]) == 1); % Make sure same labels
labels = labelTPTable;
clear labelTPTable labelPairTable

%% Count total number of objects with labelPairTable
n_obj_pair = 0;
for r  = 1:height(pixPairTable)
    for c = 1:width(pixPairTable)
        % get the table cell
        detection = pixPairTable{r, c}{1};
        [cell_r, cell_c] = size(detection);
        if cell_c == 1
            continue
        end
        n_obj_pair = n_obj_pair + size(detection{1},1);
    end
end
fprintf("Total number of objects given by pixPairTable: %d\n",n_obj_pair);

%% Count total number of objects with TP table
% There are invalid objects in labelTPTable needs to be removed.
% Class/label information is discarded.
n_obj_TP = 0;   % total number of objects
n_invalid = 0;
for r = 1:height(rectTPTableClean)
    for c = 1:width(rectTPTableClean)
        objs = rectTPTableClean{r,c}{1};   % empty array or cell array
        if size(objs,1) ~= 0
            % Find empty cells, there should be 10 of them.
            emptyCells = [];
            for i = 1:size(objs,1)
                obj = objs(i,:);      % 1 * 5 cell
                if size(obj{1}, 1) == 0
                    emptyCells = [emptyCells, i];
                end
            end
            
            % delete empty cells if found.
            if isempty(emptyCells) == 0
                rectTPTableClean{r,c}{1}(emptyCells, :) = [];
                
                % how many rows based on pixPairTable
                s1 = size(pixPairTable{r,c}{1}{1},1);
                % how many rows based on rectTPTable
                s2 = size(rectTPTableClean{r,c}{1},1);
                if s1 ~= s2
                    fprintf("object count un-match: row-%d, col-%d, ",r,c);
                    fprintf("diff:%d\n", s1-s2);
                end
                
                objs = rectTPTableClean{r,c}{1};  % re-assign
                n_invalid = n_invalid + length(emptyCells);
            end
            
            n_obj_TP = n_obj_TP + size(objs,1);
        end
    end
end

fprintf("Removed %d from rectTPTable, as they are empty.\n", n_invalid);
fprintf("Total valid objects given by rectTPTable: %d\n", n_obj_TP);

