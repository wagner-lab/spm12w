% spm12w_userdefaults
%
% spm12w userdefaults contains parameters for overriding spm12w_defaults.
% See spm12w_defaults.m for more information. 
%
% Parameters in this file will override the same parameter in spm12w_defaults.m
% and are themselves overridden by any identical parameters in the study's
% parameter files.
%
% Note that when defaults are overriden, spm12w will indicate this in the
% command window output. If you DO NOT see these messages, then spm12w is
% not overridding the defaults. The most likely reason is that the
% spm12w_userdefaults.m file (option 1 above) is not in the matlab path or,
% if using options 2 or 3, that your parameters or glm file contains a typo. 
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: November, 2014 | Updated: December, 2014
% =======1=========2=========3=========4=========5=========6=========7=========8

% spm12w user defaults 
%   uncomment the default you want to edit and change the value below.

% Begin copy of spm12w_defaults.m
% spm12w naming conventions
% def.fmri         = 'epi';            % Prefix for functional data
% def.mri          = 'anat';           % Prefix for anatomical data
% def.mristrip     = 'brain';          % Prefix for skull stirpped anatomical data
% def.mnitok       = 'native';         % Current space of data
% def.preplog      = 'preprocess.log'; % Default prepro logfile name
% def.glmlog       = 'glm_%s.log';     % Default glm logfile name. 
%                                      % Fill in glm name at creation 
% def.niceline     = [repmat('=',1,80)]; % A nice line for printing.

% spm12w logging defaults
% Logging level: 
%         1-All messages shown ([INFO],[WARNING],[EXCEPTION],[DEBUG])
%         2-Show only [INFO] & [WARNING] & [EXCEPTION]
%         3-Show only [INFO]
%         4-No messages (spm12w will run silently)
% def.loglevel     = 1;

% spm12w paths
% p.root:       The root dir for your study
% p.anadir:     The dir where user specific analyses are stored
% p.auxdir:     The dir where behavioral and other auxil data are stored.
% p.notedir:    The dir study related notes are stored
% p.prepdir:    The dir where the current preproc pipeline is stored
% p.qadir:      The dir where preprocessing qa output is stored
% p.rawdir:     The dir containing the "cleaned" raw nifti files
% p.scriptsdir: The dir where all user specific scripts are stored
% p.datadir:    The dir where preprocessing output is stored
% p.onsetsdir:  The dir where GLM onsets are stored
% p.resdir:     The dir where GLM output will be saved
% p.glmdir:     The dir name for the specifed GLM model
% p.rfxroot:    The dir root name for all 2nd level univariate tests
% p.rfxdir:     The dir name for the specified 2nd level test
% def.root       = p_.study_dir;
% def.anadir     = fullfile(def.root,'analysis',p_.username);
% def.auxdir     = fullfile(def.root, 'auxil');
% def.notesdir   = fullfile(def.root, 'notes');
% def.prepdir    = fullfile(def.root, 'prep', p_.prep_name);
% def.qadir      = fullfile(def.root, 'qa', p_.prep_name);
% def.rawdir     = fullfile(def.root, 'raw', p_.sid);
% def.scriptsdir = fullfile(def.root, 'scripts', p_.username);
% def.datadir    = fullfile(def.prepdir, p_.sid);
% def.onsetsdir  = fullfile(def.auxdir, 'onsets', p_.onsets_dir);
% def.resdir     = fullfile(def.anadir, 'results');
% def.glmdir     = fullfile(def.resdir, p_.sid, p_.glm_dir);
% def.rfxroot    = fullfile(def.resdir, 'rfx');
% def.rfxdir     = fullfile(def.rfxroot, p_.rfx_dir);
                    
% spm12w defaults for preprocessing
% Preprocessing defaults
% def.trimvols = 0;         % Trim vols off front & back of session (default=0)
%                           % can be vector [frontVols, backVols] i.e., [5,15].
%                           % or matrix (sess X vols) [5,15; 2,12].                         
% def.shuffle    = 1;         % Shuffle check
% def.despike    = 0;         % Despiking
% def.slicetime  = 1;         % Slicetime correction
% def.realign    = 1;         % Motion correction
% def.snr        = 1;         % Calculate snr
% def.slices     = 1;         % Calculate slice noise
% def.matkill    = 1;         % Delete intermediate mat files during cleanup
                            
% Despiking defaults
% def.cut        = [2.5,4];   % Afni defaults for despiking
% def.c_order    = 'nvols/30';% Afni defaults for despiking
% def.polyord    = 2;         % Number of polynomials for despiking
% def.least_sq   = 0;         % If 0 then do L1 regression for despiking
% def.dilate     = 4;         % Dilates spm's "automask" by factor during despiking                          
                       
% Slicetime correction defaults
% def.sformula   = 'philips'; % Slice order formula ('philips' or 'regular')
% def.refslice   = 1;         % Reference slice for slice timing

% Realign defaults
% def.realign_rtm = 1;      % Realign to 1st only (0) or 2pass(1st then sess mean) (1)
% def.realign_quality = 1;  % Quality versus speed trade-off. SPM12 default 0.9

% Cleanup defaults
% def.cleanup    = 1;         % 0 = keep all files (useful for debugging)
%                             % 1 = delete all but the smoothed epi
%                             % 2 = delete all but the smoothed and unsmoothed epi
% def.cleanupzip = 0;         % gzip final datasets (0 = no). 

% Dartel defaults
% def.biasreg       = 0.0001;   % Intensity bias regularisation
% 							  % (0 | 0.00001 | *0.0001 | 0.001 | 0.01 | 0.1 | 1 | 10) 
% def.mrf           = 2;        % Segmentation cleanup via MRF (0 | 1 | *2 | 3 | 4)
% def.warpreg       = 4;        % Warp regularisation (spm default:4|can use 0.4 for T1)                        
% def.sampling      = 3;        % Sampling distance (1|2|*3) < 3mm will slow computation 
% def.tissues       = 2;        % Tissue classes for template (1:GM|*2:GM+WM|3:GM+WM+CSF)
% def.dartelsmooth  = 6;        % Smooting kernel for DARTEL to MNI (GM,WM,EPI only)

% Normalization and write parameters
% Interpolation set to 7 bsplines gives a less blury normalized EPI.
% After smoothing that extra spatial resolution is lost.
% However, for studies with low smoothing, it might be worth bumping up
% the interp_r and interp_w and interp_c.
% r = realign, w = unwarp, c = coregister
% def.wrap_r     = [0 0 0];     %SPM12/8 default [0 0 0] | spm2w was [0 1 0]
% def.wrap_w     = [0 0 0];     %SPM12/8 default [0 0 0] | spm2w was [0 1 0]
% def.interp_r   = 7;           %SPM12/8 default 4       | spm2w was 7
% def.interp_w   = 7;           %SPM12/8 default 4       | spm2w was 7
% def.interp_c   = 7;           %SPM12/8 default 4       | spm2w was 4
% def.voxsize    = [3 3 3];     %SPM12/8 default [2x2x2] | Don't resample
% def.boundbox   = [-78 -117 -72;...
%                  78 76 86];   %SPM12 default is [-78 -112 -70; 78 76 85]
%                               %SPM8 default is [-78 -112 -50; 78 76 85] 
% 							  %But that's for 2x2x2.
%                               %Ours works better for 3x3x3
                              
% Brain mask
% Can use either bigmask in 1x1x1 or 3x3x3
% def.mask = 'bigmask.nii';