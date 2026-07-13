%% PorousLab - FEM framework for multiphysics problems in porous media
%
% This is the regression testing file of PorousLab.
% To run tests, execute this file to compare current results with reference values
% considering a given tolerance.
%
clc; clearvars; close all;
tol = 1e-8;

files_dir = fullfile(pwd, 'files');
test_files = dir(fullfile(files_dir, 'test_*.m'));
src_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src');
addpath(genpath(src_dir));

total = 0;
passed = 0;
failed = 0;
skipped = 0;

for i = 1:length(test_files)
    test_name = test_files(i).name(1:end-2);
    fprintf('Running %s...\n', test_name);

    clear mdl;
    run(fullfile(files_dir, test_name));
    total = total + 1;

    ref_file = fullfile(files_dir, [test_name,'_ref.mat']);
    if ~exist(ref_file, 'file')
        warning('⚠️ Reference file missing: %s\n\n', ref_file);
        skipped = skipped + 1;
        continue;
    end

    load(ref_file, 'U_ref');

    if ~isequal(size(mdl.U), size(U_ref))
        fprintf('❌ Failed: size mismatch (mdl.U is [%s], U_ref is [%s])\n\n', num2str(size(mdl.U)), num2str(size(U_ref)));
        failed = failed + 1;
        continue;
    end

    err = norm(mdl.U - U_ref) / norm(U_ref);

    if err < tol
        fprintf('✅ Passed (rel. error = %.2e)\n\n', err);
        passed = passed + 1;
    else
        fprintf('❌ Failed (rel. error = %.2e)\n\n', err);
        failed = failed + 1;
    end
end

fprintf('======== Test Summary ========\n');
fprintf('Total tests:........%d\n', total);
fprintf('Passed:.............%d\n', passed);
fprintf('Failed:.............%d\n', failed);
fprintf('Skipped (no ref):...%d\n', skipped);
fprintf('===============================\n');
