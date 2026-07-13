%% obtainHTML Function
% This script generates HTML documentation for all MATLAB scripts (.m files)
% within the porousLab repository, excluding specific folders and files.
% The generated HTML files are saved in the _docs/html_ directory of the repository.
% If _docs/html_ doesn't exist, it is created.
% 
%% Functionalities
% 1. Identifies the current script's path and repository structure.
% 2. Recursively searches for all .m files in the repository.
% 3. Excludes specified folders and files from processing.
% 4. Publishes each script to HTML format using MATLAB's _publish_ 
%    function.
% 5. Saves the generated HTML files in the _docs/html_ directory.
% 6. Handles errors during publishing and logs any failed scripts.
% 
%% Usage
% Run this script from within the porousLab repository to generate HTML 
% documentation for all MATLAB scripts, excluding specified folders and files.
% 
%% Specify input and output paths
clear; clc; close all;

% Get the path and name of the current script
currentScriptPath = mfilename('fullpath');
currentScriptName = strcat(mfilename, '.m'); % Add .m extension
currentFolder = fileparts(currentScriptPath);

% Move up two levels to reach the root of the repository (porousLab)
repoFolder = fileparts(fileparts(currentFolder));

% Define the output folder inside porousLab/docs/html
outputFolder = fullfile(repoFolder, 'docs', 'html'); 

% Ensure output folder exists
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Obtain all files in the repository following the specified criteria
% Get a list of all .m files in the porousLab repository (all subdirectories)
mFiles = dir(fullfile(repoFolder, '**', '*.m')); % Recursively finds all .m files

% Exclude specific directories and files
excludedFolders = {'examples','tests'}; % Add more folders if needed
excludedFiles = {'print_header'}; % Add more filenames if needed

% Create the filtered files list
filteredFiles = [];

% List to keep track of scripts that failed
failedScripts = {};

% Filter the files considering the excludedFolders and the excludedFiles
for k = 1:length(mFiles)
    relativePath = erase(mFiles(k).folder, repoFolder); % Get relative path
    isExcludedFolder = any(contains(relativePath, excludedFolders)); % Check exclusion for folder
    isExcludedFile = any(strcmp(mFiles(k).name, excludedFiles)); % Check exclusion for file
    isCurrentScript = strcmp(mFiles(k).name, currentScriptName); % Check if it is the current script
    
    if ~isExcludedFolder && ~isExcludedFile && ~isCurrentScript
        filteredFiles = [filteredFiles; mFiles(k)]; %#ok<AGROW>
    end
end

%%  Filter any repeated file name
% Obtain the names
fileNames = {filteredFiles.name};

% Find unique file names and their indices
[uniqueFileNames, uniqueIndices] = unique(fileNames);

% Find duplicate file names by comparing the indices
duplicateIndices = setdiff(1:length(fileNames), uniqueIndices);

% If there are duplicates, display them
if ~isempty(duplicateIndices)
    fprintf('====== Duplicated file names ======\n');
    for idx = 1:length(duplicateIndices)
        disp(['    ' fileNames{duplicateIndices(idx)}]);
        disp(' ')
    end
else
    disp('No duplicate file names found.');
    disp(' ')
end

%% Publish the files
% Filter the filteredFiles list to include only non-duplicate files
nonDuplicateFiles = filteredFiles(uniqueIndices);

% Start generating HTML files message
fprintf('===== HTML generation started =====\n')

% Publish each non-duplicate file to the specified output folder
for k = 1:length(nonDuplicateFiles)
    scriptPath = fullfile(nonDuplicateFiles(k).folder, nonDuplicateFiles(k).name);
    
    % Temporarily add the script's folder to the MATLAB path
    addpath(nonDuplicateFiles(k).folder);
    
    options = struct('format', 'html', 'outputDir', outputFolder, 'evalCode', false);
    try
        publish(scriptPath, options);
        fprintf('Published (%d/%d): %s\n', k, length(nonDuplicateFiles), nonDuplicateFiles(k).name);
    catch ME
        fprintf('Failed to publish: %s\nError: %s\n', nonDuplicateFiles(k).name, ME.message);
        % Add failed script to the list
        failedScripts{end+1} = nonDuplicateFiles(k).name;
    end
    
    % Remove the script's folder from the MATLAB path to avoid conflicts
    rmpath(nonDuplicateFiles(k).folder);
end

% End generating HTML files message
fprintf('====== HTML generation ended ======\n')

%% Final info
% Final info message
disp(' ')
fprintf('============ Final info ===========\n')

% After publishing, display or save the list of failed scripts
if ~isempty(failedScripts)
    fprintf('The following scripts failed to publish:\n');
    disp(failedScripts);
else
    fprintf('All scripts were published successfully!\n');
end

% Final output message
fprintf('(%d/%d) scripts have been published to: %s\n', k, length(filteredFiles), outputFolder);
