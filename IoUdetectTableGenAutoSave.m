%% This creates a table of all detected images over the a given dataset
clear

%% Get all csv filenames from the desired directory
Files=dir('..\All_Combined\Detected_CSV_Results\*.csv');
filesSorted = natsortfiles({Files.name});   % sort files in chronological order
saveDir = '.\results';       % Current folder
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

%% get class labels. These should be consistent across all data
load(".\ThesisLabels.mat");   % labelDefs, comes from the previous section
classes = labelDefs.Name';  % 1 x 10 cell, 10 classes
classes = flip(classes);    % may not necessary. Be consistent with Ben's original setting.

% For some unknown reason, if the table is not initialized with an empty row of cells to start, 
% the program will crash if more than one object is assigned to a class.
% this row is then deleted at the end of the program. This is really dumb.
IoUdetected = cell2table(cell(1,length(classes)+1),'VariableNames', ['imageFilename', classes]);
IoUdetectedWithConfidence = cell2table(cell(1,length(classes)+1),'VariableNames', ['imageFilename', classes]);

%% Read each file
fprintf("Processing %d result csv files:", length(Files));
for k=1:length(Files)
    % print progress info
    if rem(k, 100) == 1
        fprintf("%d, ", k);
    end
    [IoUdetectedInit, IoUdetectedWithConfidenceInit] = openCSV_fixed(filesSorted{k}, Files(k).folder, classes);
    IoUdetected = [IoUdetected; IoUdetectedInit];
	IoUdetectedWithConfidence = [IoUdetectedWithConfidence; IoUdetectedWithConfidenceInit];
end
fprintf("done.\n")

%% Save results in one file
IoUdetected(1,:)=[];
IoUdetectedWithConfidence(1,:)=[];
save([saveDir,'\IoUdetected'],'IoUdetected')
save([saveDir,'\IoUdetectedWithConfidence'],'IoUdetectedWithConfidence')
fprintf("Results saved in %s: \n\tIoUdetected.mat \n\tIoUdetectedWithConfidence.mat\n\n", saveDir)
