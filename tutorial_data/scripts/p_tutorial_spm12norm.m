% spm12w r6906
% Parameters file for fMRI preprocessing
% Last updated: September 2016
% =======1=========2=========3=========4=========5=========6=========7=========8

% User name
p.username = 'ddw';

% Preprocessing output directory
p.prep_name = 'spm12_norm';

% Preprocessing Routines - 1=yes 0=no
p.unwarp        = 1; % Unwarping (correct field inhomogeneties)      
p.normalize     = 'spm12'; % Normalize type ('none','epi','spm12','dartel')
p.smoothing     = 8; % Size of smoothing kernel in FWHM (0 for no smoothing)