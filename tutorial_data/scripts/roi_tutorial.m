% spm12w r6225
% Parameters file for roi analysis
% Last updated: May, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% User name
roi.username = 'ddw';

% Paths and names
roi.study_dir = 'C:\Users\ddw\Documents\MATLAB\tutorial_data';
roi.glm_name = 'glm_tutorial';      % glm to use for extracting parameters
roi.roi_name = 'roi_tutorial';      % name of directory for roi analysis output

% CSV files specifying roi specs and subject variables (optional)
roi.spec_file = 'tutorial_spec.csv'; % name of a csv file specifying rois
roi.var_file = 'tutorial_vars.csv';  % name of a csv file with subject variables

% 1st level GLM contrasts to use for roi analysis 
roi.conds = {'allVSbaseline','humVSall'};

% ROI specifications - (X,Y,Z, sphere diameter in mm) or mask filename
%                    - for masks use nifti image name instead of coordinate
%                    - if no sphere diameter specified, default will be used
%                    - can be used in conjunction with an roi spec file. 
roi.roi.l_dmpfc_ba10 = [-3,63,15,6];
roi.roi.l_precuneus_ba32 = [-3,-60,24];
roi.roi.l_amygdala = 'aal_l_amyg_3x3x3.nii';

% ROI statistics - statistics: descriptives, ttest1 ttest2, correl1, correl2
%                - leave blank or omit if statistic not desired
%                - use cell arrays for stats on multiple contrasts
%                - use strings to define a contrast formula
%                - 'all_conditions' is a reserved word
roi.stats.descriptives = 'all_conditions';
roi.stats.ttest1 = {'allVSbaseline','humVSall vs. allVSbaseline'};
roi.stats.ttest2 = {'allVSbaseline', 'allVSbaseline humVSall'};
roi.stats.correl1 = 'humVSall';
roi.stats.correl2 = 'allVSbaseline';