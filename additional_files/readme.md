spm12w
======
The files in this directory should be added to your locally installed copy of spm12.

canonical.7z
------------
The archive canonical.7z contains high quality brains in MNI space that may be used
for visualization. These files need to be placed inside the spm12/canonical directory.

spm12
-------
Dead simple bash script to start matlab in terminal mode and run the spm12 path script.
Only useful for Linux and possibly Mac.

spm12path.m
----------------
Example script of how to setup your matlab paths for spm12 and spm12w. Although you can
manually add spm12/spm12w and its associated directories to your Matlab path via the gui,
this method is preferable as it allows multiple versions of spm12 and multiple branches
of spm12w to coexist (simply change the paths in this script to point to other versions).
