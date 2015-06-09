% spm12w r6225
% Parameters file for fMRI preprocessing
% Last updated: October, 2014
% =======1=========2=========3=========4=========5=========6=========7=========8

% User name
p.username = 'ddw';

% Paths and names
p.study_dir = 'C:\Users\ddw\Documents\MATLAB\tutorial_data';
p.prep_name = 'standard_preprocessing';
p.rfx_name  = 'standard_preprocessing';

% Preprocessing Routines - 1=yes 0=no
p.unwarp        = 1; % Unwarping (correct field inhomogeneties)      
p.normalize     = 'epi'; % Normalize type ('epi','spm12','dartel')
p.smoothing     = 8; % Size of smoothing kernel in FWHM (0 for no smoothing)

% Segmentation, Dartel and VBM8 parameters - 1=yes 0=no
p.coreg2epi     = 1;        % Coregister Anat to EPI data (mean of session 1)?
p.whichseg      = 'DARTEL'; % Segmentation type? (DARTEL | VBM8)
p.whichtemplate = 'DARTEL'; % Template style? VBM8 uses IXI550 (DARTEL | VBM8)
p.whichnorm     = 'DARTEL'; % Target for EPI normalization? (DARTEL | VBM8)
