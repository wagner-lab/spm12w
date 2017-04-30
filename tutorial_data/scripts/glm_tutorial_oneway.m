% spm12w r6906
% Parameters file for 1st level glm analysis and 2nd level rfx analysis
% Last updated: March, 2017
% =======1=========2=========3=========4=========5=========6=========7=========8
%
% Tutorial note: 
% The purpose of this file is to give you an exmaple of how to run a simple
% GLM and RFX analysis and how to specify the various types of contrasts 
% for a simple event-related study.
%
% Three ways to run this file: 
% 1> spm12w    #Use the gui and select "Estimate GLM"
% 2> spm12w_glm_compute('sid','s01','glm_file','glm_tutorial.m')
% 3> spm12w('stage','glm', 'sids','s01','para_file','glm_tutorial.m')
%
% After estimation, you can inspect the model in SPM (SPM -> REVIEW) or view
% the PDF file in the GLM output directory. The next step is to generate
% contrasts for each subject followed by a random effects analysis.
%
% Final note: Additional parameters are available in spm12w_defaults.m. For
% example, the glm explicit mask, the microtime resolution and onsets, switches
% to exclude motion regressors, use a high-pass-filter etc. To override a
% default, copy it to this file and change the structure name from def.<field>
% to glm.<field>. For an example of a default override, see the
% glm_tutorial_hpf_ar1.m file ovverrides the hpf and autocorrelation correction
% defaults. 

% User name
glm.username = 'ddw';

% GLM input/output directory
glm.prep_name = 'epi_norm';
glm.glm_name = 'glm_tutorial';

% GLM onsets directory name
glm.ons_dir  = 'tutorial';

% GLM Conditions (cell arrays seperated by commas)
glm.events     = {'human','animal','vegetable','mineral'};
glm.blocks     = {};
glm.regressors = {};

% GLM Contrasts - Numeric contrasts should sum to zero unless vs. baseline
%               - String contrasts must match Condition names.
%               - String contrast direction determine by the placement of 'vs.'
%               - Special keywords: 'housewine' | 'fir_bins' | 'fir_hrf'
% Note that the housewine's allVSbaseline will skip pmods but is not savy
% to complex designs (i.e., where you want all conditions VS baseline but
% not your user supplied regressors under glm.regressors, so in those cases
% you should create your own allVSbaseline and refuse the housewine). 
glm.con.housewine = 'housewine';
glm.con.humVSveg  = [1 0 -1 0];
glm.con.aniVSmin  = [0 1 0 -1];
glm.con.minANDhum = 'mineral human';
glm.con.humVSall  = 'human vs. animal vegetable mineral';

% RFX Specification
% Note that we "should" put the rfx specifications in a seperate file from the
% glm, however to cut down on the number of parameter files and for historical
% reasons (e.g., spm2w, spm8w) we tack the rfx specifications here.

% The assignemtns below will run a single one-way repeated measures ANOVA
% across the conditions humVSbaseline, animalVSbaseline, vegetableVSbaseline
% and mineralVSbaseline. The contrasts were not created by the user, but are
% part of the free drinks you get when you ask for the housewine contrast.
% Here, the assignments opperate differently than with the one-sample test 
% used in the other glm tutorial files. For instance, in the one-sample case,
% the glm.rfx_conds specifies on which conditions the one-sample t-test will
% be performed and the output of these analysis will be written to directories
% with the same name as the conditions with the directory name specified in
% glm.rfx_name) serving as the parent directory. Here, the oneway_conds variable 
% specifies which conditions0 to include in a single one-way ANOVA. 
% The name of the directory in which the anova results will be written to is 
% given by glm.rfx_name. 

glm.rfx_type  = 'anova1'; % Normally we don't include type as we usually do 
                          % one-sample and that's the default. Here we have to 
                          % explicitly specify the type.
glm.rfx_name = 'oneway_tutorial';
glm.rfx_conds = {'humanVSbaseline', 'animalVSbaseline', ...
                 'vegetableVSbaseline', 'mineralVSbaseline'};