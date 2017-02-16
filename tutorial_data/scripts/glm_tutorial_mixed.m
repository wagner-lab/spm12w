% spm12w r6225
% Parameters file for 1st level glm analysis and 2nd level rfx analysis
% Last updated: October, 2016
% =======1=========2=========3=========4=========5=========6=========7=========8
%
% Tutorial note: 
% The purpose of this file is to give you an exmaple of how to run a mixed
% block-event design (i.e., State-Item design) following the Wash U method.
% Specifically: To avoid state's stealing item effets, blocks are not
% convolved with the HRF. Item effects are specified usign FIR. 
%
% Three ways to run this file: 
% 1> spm12w    #Use the gui and select "Estimate GLM"
% 2> spm12w_glm_compute('sid','s01','glm_file','glm_tutorial.m')
% 3> spm12w('stage','glm', 'sids','s01','para_file','glm_tutorial.m')

% User name
glm.username = 'ddw';

% GLM input/output directory
glm.prep_name = 'epi_norm';
glm.glm_name = 'glm_tutorial_mixed';

% GLM onsets directory name
glm.ons_dir  = 'tutorial_mixed';

% GLM Conditions (cell arrays seperated by commas)
glm.events     = {'human','animal','oncue','offcue'};
glm.blocks     = {'implicit','explicit'};
glm.regressors = {'psychophys'};

% GLM override defaults
glm.block_conv = 0;     % Convolve blocks (0 for state/item) - 1=yes 0=no

% GLM HRF Specifications:
% OPTIONS:'hrf'
%         'hrf (with time derivative)'
%         'hrf (with time and dispersion derivatives)'
%         'Fourier set'
%         'Fourier set (Hanning)'
%         'Gamma functions'
%         'Finite Impulse Response'
glm.hrf       = 'Finite Impulse Response';  %Set to Finite Impulse Response for MIXED designs.
glm.hrfwindow = 20;     %For FIR: Length of HRF window 
glm.hrfbasis  = 8;      %For FIR: Number of bins per window

% GLM Contrasts - Numeric contrasts should sum to zero unless vs. baseline
%               - String contrasts must match Condition names.
%               - String contrast direction determine by the placement of 'vs.'
%               - Special keywords: 'housewine' | 'fir_bins' | 'fir_hrf'
% Note that the housewine's allVSbaseline will skip pmods but is not savy
% to complex designs (i.e., where you want all conditions VS baseline but
% not your user supplied regressors under glm.regressors, so in those cases
% you should create your own allVSbaseline and refuse the housewine). 
glm.con.housewine = 'housewine';
glm.con.socVSnsoc = 'implicit vs. explicit';
glm.con.psychoreg = 'psychophys';

% RFX Specification
glm.rfx_name = 'rfx_tutorial';
glm.rfx_conds = {'allVSbaseline','humVSall','socVSnsoc','psychoreg'};