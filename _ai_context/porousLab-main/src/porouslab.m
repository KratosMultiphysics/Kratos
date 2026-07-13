%% PorousLab - FEM framework for multiphysics problems in porous media
%
% This is the main function of PorousLab.
% To run a simulation, call this function by passing the full path to the simulation script as input argument.
% If no input argument is provided, a dialog appears for manual selection of the simulation script.
%
function mdl = porouslab(varargin)
    clc; clearvars -except varargin; close all;
    addpath(genpath(pwd));
    print_header;

    if nargin == 0
        [file_name, file_path] = uigetfile('*.m', 'Select a script to run');
        if isequal(file_name, 0)
            disp('No file selected.');
            return;
        end
    elseif nargin == 1
        script = varargin{1};
        if ischar(script) || isstring(script)
            script = char(script);
            [file_path, file_name, ext] = fileparts(script);
            file_name = [file_name, ext];
            if ~strcmp(ext,'.m') || ~exist(script,'file')
                error('The provided file must be a valid Matlab script.');
            end
        else
            error('Input must be a string or character vector.');
        end
    else
        error('Invalid number of input arguments.');
    end

    addpath(file_path);
    run(file_name);

    % Get output variables
    vars = whos;
    ws = struct();
    for k = 1:length(vars)
        varName = vars(k).name;
        ws.(varName) = eval(varName);
    end
end
