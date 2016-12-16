% spm12w r6225
% Parameters file for 1st level glm analysis and 2nd level rfx analysis
% Last updated: October, 2016
% =======1=========2=========3=========4=========5=========6=========7=========8

% Tutorial note: 
% The purpose of this file is to give you an example of how to generate a GLM
% model without any data. This is useful in instances where you are designing
% a study and would like to see what the GLM model will look in order to inspect
% inpsect it for collinearity between conditions.
%
% Two ways to run this file: 
% 1> spm12w_glm_compute('sid','s01','glm_file','glm_tutorial_design_only.m')
% 2> spm12w('stage','glm', 'sids','s01','para_file','glm_tutorial_design_only.m')

% User name
glm.username = 'ddw';

% GLM output directory
glm.prep_name = 'epi_norm';
glm.glm_name = 'glm_tutorial_designonly';

% GLM onsets directory name
glm.ons_dir  = 'tutorial';

% GLM Onsets File Specifications
glm.design_only = 1;   % Design only (i.e., no data)
glm.outliers    = 0;   % Include outlier as nuissance in GLM  

% Tutorial Note: 
% Because we are estimate the GLM only nad have no data, you must manually
% specify the intended number of sessions, nvols, TR and disable including
% movement regressors.

glm.nses = 2;
glm.nvols = [120,120];
glm.tr = [2.5,2.5];
glm.move = 0;
                       
% GLM Conditions (cell arrays seperated by commas)
glm.events     = {'human','animal','vegetable','mineral'};
glm.blocks     = {};
glm.regressors = {};