% spm12w r6906
% Parameters file for 1st level glm analysis and 2nd level rfx analysis
% Last updated: March, 2017
% =======1=========2=========3=========4=========5=========6=========7=========8
%
% Tutorial note: 
% The purpose of this file is to give you an exmaple of how to run a simple
% GLM and RFX analysis and how to specify the various types of contrasts 
% for a simple event-related study with a custom microtime resolution and onset.
%
% Three ways to run this file: 
% 1> spm12w    #Use the gui and select "Estimate GLM"
% 2> spm12w_glm_compute('sid','s01','glm_file','glm_tutorial_microtime.m')
% 3> spm12w('stage','glm', 'sids','s01','para_file','glm_tutorial_microtime.m')
%
% After estimation, you can inspect the model in SPM (SPM -> REVIEW) or view
% the PDF file in the GLM output directory. The next step is to generate
% contrasts for each subject followed by a random effects analysis.
%
% Here we override the default microtime resolution(16) and onset(1). This is
% useful in situation where your slice time correction reference slice does
% *not* correspond to the first acquisition slice. For instance, if your slice
% acquisition order starts with slice 7 but you decide to use slice 1 as the
% slicetime correction reference slice, you will want to set your microtime
% onset to the timebin that is associated with slice 1. Note that during
% slice time correction, the preprocessing script will recommend a 
% microtime resolution and onset for you to use durin GLM estimation. The 
% preprocessing log file will contain these values.

% User name
glm.username = 'ddw';

% GLM input/output directory
glm.prep_name = 'epi_norm';
glm.glm_name = 'glm_tutorial_microtime';

% GLM onsets directory name
glm.ons_dir  = 'tutorial';

% GLM Conditions (cell arrays seperated by commas)
glm.events     = {'human','animal','vegetable','mineral'};
glm.blocks     = {};
glm.regressors = {};

% GLM Time speifications: 
glm.tbins     = 36; % Time bin resolution for each scan (16, SPM default)
glm.tref      = 1;  % Reference time bin for GLM. Usually timebin corresponding 
                    % to the ref slice during slicetime correction. For Philips
                    % interleaved the first slice acquired is slice 1 so the
                    % the reference timebin is also 1 if p.refslice=1. 

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
glm.rfx_name = 'rfx_tutorial_microtime';
glm.rfx_conds = {'allVSbaseline','humVSall'};