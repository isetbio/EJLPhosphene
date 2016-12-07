function rootPath=phospheneRootPath()
% Return the path to the root phosphene directory
%
% This function must reside in the directory at the base of the phosphene
% directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Example:
%   fullfile(isetbioRootPath,'data')

rootPath=which('phospheneRootPath');

[rootPath,~,~]=fileparts(rootPath);

return
