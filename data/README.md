# data/

Nothing in this folder is tracked by git (see `.gitignore`).

Expected input: either a `workspace.mat` written by
`hsi-chromatic-defocus-depth/src/sparse/MainScript.m` (holds the normalised
cube `normat` and the sparse depth), or SnapScan cubes saved as `.mat` files.
The scenes used in the thesis are lab set-ups with objects partly hidden
behind others.

The original scripts had absolute paths to the author's machines. These were
replaced with relative placeholders that start with `.\data\...`, so the
simplest set-up is to copy (or symlink) your data here, keeping the same
sub-folder names, and run scripts from the repository root. You can also just
edit the path variables at the top of each script.
