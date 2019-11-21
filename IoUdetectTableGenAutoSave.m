% This creates a table of all detected images over the a given dataset

clear
% Get all csv filenames from the desired directory

Files=dir('C:\Users\benmi\Documents\Thesis\Thesis Instructions\Raw Data\Detected_CSV_Results\*.csv');
saveDir = 'C:\Users\benmi\Documents\Thesis\Thesis Instructions\Raw Data';

% Don't change anything below here!!!

% sort files into chronological order
filesSorted = natsortfiles({Files.name});
% get class labels. These should be consistent across all data
load('C:\Users\benmi\Documents\Thesis\Matlab\IoU - Intersection Over Union and Performance Metrics\groundTruthFrame1Test2.mat');
% load("C:\Users\benmi\Documents\Thesis\Training Results\SnowAndNight 4 Class Results\testLabels.mat");
classes = {IoUtruth.Properties.VariableNames{:}};
% classes = {testLabelData.Properties.VariableNames{:}};
classes = classes(2:end);
% IoUdetected = cell2table(cell(0,length(classes)+1),'VariableNames', ['imageFilename', classes]);

%For some unkown reason, if the table is not initialized with an empty row of cells to start, 
% the program will crash if more than one object is assigned to a class.
% this row is then deleted at the end of the program. This is really dumb.
IoUdetected = cell2table(cell(1,length(classes)+1),'VariableNames', ['imageFilename', classes]);

% IoUdetectedWithConfidence = cell2table(cell(0,length(classes)+1),'VariableNames', ['imageFilename', classes]);
IoUdetectedWithConfidence = cell2table(cell(1,length(classes)+1),'VariableNames', ['imageFilename', classes]);


for k=1:length(Files)
%    [IoUdetectedInit, IoUdetectedWithConfidenceInit] = openCSV(Files(k).name, Files(k).folder, classes);
   [IoUdetectedInit, IoUdetectedWithConfidenceInit] = openCSV_fixed(filesSorted{k}, Files(k).folder, classes);
   IoUdetected = [IoUdetected; IoUdetectedInit];
   IoUdetectedWithConfidence = [IoUdetectedWithConfidence; IoUdetectedWithConfidenceInit];
end

IoUdetected(1,:)=[];
IoUdetectedWithConfidence(1,:)=[];
save([saveDir,'\IoUdetected'],'IoUdetected')
save([saveDir,'\IoUdetectedWithConfidence'],'IoUdetectedWithConfidence')


