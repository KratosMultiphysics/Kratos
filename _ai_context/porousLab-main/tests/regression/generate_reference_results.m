%% PorousLab - FEM framework for multiphysics problems in porous media
%
% This is the reference generation file for regression tests in PorousLab.
% Execute this file to generate or overwrite the reference results of each test case.
%
clc; clearvars; close all;

files_dir = fullfile(pwd, 'files');
test_files = dir(fullfile(files_dir, 'test_*.m'));
src_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src');
addpath(genpath(src_dir));

for i = 1:length(test_files)
    test_name = test_files(i).name(1:end-2);
    fprintf('Running %s...\n', test_name);

    try
        clear mdl U_ref;
        run(fullfile(files_dir, test_name));

        if ~exist('mdl', 'var') || ~isprop(mdl, 'U') || isempty(mdl.U)
            warning('⚠️ Skipped: mdl.U not found\n\n');
            continue;
        end

        U_ref = mdl.U;
        ref_name = [test_name,'_ref.mat'];
        save(fullfile(files_dir, ref_name), 'U_ref');
        fprintf('✅ Saved reference results in %s\n\n', ref_name);

    catch ME
        warning('⚠️ Skipped: Failed to run\n\n');
        continue;
    end
end
