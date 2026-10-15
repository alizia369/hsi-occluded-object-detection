function setup_paths()
%SETUP_PATHS Add this repository's code folders to the MATLAB path.
%   Run once per MATLAB session from anywhere:  >> setup_paths
%   Data and result folders are left off the path on purpose.

root  = fileparts(mfilename('fullpath'));
parts = strsplit(genpath(root), pathsep);
skip  = {[filesep '.git'], [filesep 'data'], [filesep 'results'], [filesep 'Results']};
keep  = true(size(parts));
for k = 1:numel(skip)
    keep = keep & cellfun(@isempty, strfind(parts, skip{k}));
end
parts = parts(keep & ~cellfun(@isempty, parts));
addpath(strjoin(parts, pathsep));
fprintf('Added %d folders from %s to the MATLAB path.\n', numel(parts), root);
end
