function output_table = runProfile
% Function to run read the profile data from the dataray camera and then
% automatically run the beam fitting routine. 

%% Profile specific parameters
% Update the wavelength and translation distance to match the current
% profile specifications

% translation distance (for shifting distances to put zero on atoms)
translate_dist = -12.4;

% specify wavelength
lambda = 532e-7;

%% Standard functionality
% This is the routine to read in files then run BeamProfileFits. This part
% of the process shouldn't change day to day

% Add subfolder 'Library' to path for internal functions
addpath(genpath([pwd filesep 'Library']));

% Get filenames of data in directory
filenames = struct2table(dir(['dataray' filesep '*.log']));

% Intialize loop
output_table = table;
for i = 1:height(filenames)
    % Get the distance from the filename
    dist = strsplit(filenames.name{i},'.log');
    dist = str2double(dist{1}) + translate_dist;
    
    % load in data and put data into the table
    output_table = fill_row(output_table, dist, importDataray([pwd filesep 'dataray' filesep filenames.name{i}]));
end

BeamProfileFits_v5_Weighted('rawData',output_table(:,1:5), 'lambda', lambda)

% Restore to original path
rmpath([pwd filesep 'Library' filesep 'Archive']);
end