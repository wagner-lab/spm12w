% spm12w r6225
% Parameters file for roi analysis
% Last updated: April, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% User name
roi.username = 'ddw';

% Paths and names
roi.study_dir = 'C:\Users\ddw\Documents\MATLAB\tutorial_data';
roi.glm_name = 'glm_tutorial';      % glm to use for extracting parameters
roi.roi_name = 'roi_tutorial';      % name of directory for roi analysis output

% Files specifying roi specifications and subject variables (optional)
roi.spec_file = 'roi_tutorial.csv'; % name of a csv file specifying rois
roi.var_file = 'roi_tutorial.csv';  % name of a csv file with subject variables

% 1st level GLM contrasts to use for roi analysis 
roi.conds = {'allVSbaseline','humVSall'};

% ROI specifications - (mni coordinates, sphere diameter in mm)
%                    - for masks use nifti image name instead of coordiante
roi.roi.dmpfc_ba10 = {[-3,63,18],6};
roi.roi.precuneus_ba32 = {[3,-60,24],6};
roi.roi.amygdala = {'aal_l_amyg_3x3x3.nii'};

% ROI statistics - statistics: descriptives, ttest1 ttest2, correl1, correl2
%                - leave blank if statistic not desired
%                - use cell arrays for stats on multiple contrasts
%                - use strings to define formula
%                - 'all_conditions' is a reserved word
roi.stats.descriptives = 'all_conditions';
roi.stats.ttest1 = {'allVSbaseline','humVSall vs. allVSbaseline'};
roi.stats.ttest2 = {'allvsbaseline, humVSall'};
roi.stats.correl1 = 'humVSall';
roi.stats.correl2 = 'allVSbaseline';