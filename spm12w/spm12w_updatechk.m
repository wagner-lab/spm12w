function spm12w_updatechk()
% spm12w_updatechk()
% 
% Simple script to print out the revision number of various core SPM
% files used by spm12w in order to faciliate determining which need
% to be combed over for modification into spm12w.
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: March, 2014 | Updated: December, 2014
% =======1=========2=========3=========4=========5=========6=========7=========

% Get paths
spm12_path = fileparts(which('spm.m'));
spm12w_path = fileparts(which('spm12w.m'));
addpath([spm12_path,'/toolbox/OldNorm/']);
addpath([spm12_path,'/toolbox/DARTEL/']);

edit_list = {'realign','slice_timing','uw_estimate','normalise',...
    'normalise_disp','preproc_run','maff8','preproc8','preproc_write8',...
    'dartel_template','dartel_norm_fun','dartel_warp'};

% Print version match
[tmp,revision] = spm('Ver');
fprintf('spm12 revision: r%s\n', revision);
for i = 1:length(edit_list)
    file12  = sprintf('spm_%s.m',edit_list{i}); 
    file12w = sprintf('spm12w_%s.m',edit_list{i}); 
    r12 = spm('Ver',file12);
    r12w = spm('Ver',file12w);
    fprintf('Revision: r%s | spm12 file: %s\n', r12, file12);
    fprintf('Revision: r%s | spm12w file: %s\n', r12w, file12w);
end