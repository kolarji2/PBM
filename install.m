%% Installation of PBM
% Step 1 - Extract all files or clone repository from github
%
% Step 2 - Usage
%
% a) Temporary usage
%       1) Copy folders @PBM, @Utils and @OptimizePBM to your
%          project and use them.

%       2) Or you can temporary add class folders to matlab path from any
%          project:
%
%               addpath(PATH_TO_PBM_FOLDER)
%
%          where in PATH_TO_PBM_FOLDER are located folders @PBM, @Utils
%          and @OptimizePBM
%
% b) Permanent usage - permanently add to matlab path
%
%   1) Make sure this script, install.m is in the same directory
%      as directories @PBM, @Utils and @OptimizePBM
%   2) Run this script
%
%   3) If you would like to remove the package go to matlab "Set Path" in
%   Environment tab. Locate link to this package. Remove and save.

if exist('@PBM','dir') && exist('@Utils','dir') && exist('@OptimizePBM','dir')
    % Temporary add script files
    addpath(pwd);    
    
    % Confirm if we should add them permanently 
    b = questdlg('Do you want to permanently add PBM class to your Matlab path?','Installation of PBM','Yes','No','Yes');
    
    if strcmp(b,'Yes')
        save = savepath;
        if save == 0
            helpdlg('PBM class was permanently added to your Matlab path.')
        else
            errordlg('The matlab path could not be saved. Check if you have permissions to do that.','Error saving matlab path');
        end
    end
    
else
    errordlg('Folders @PBM, @Utils or @OptimizePBM were not found.','Invalid directory.');
end