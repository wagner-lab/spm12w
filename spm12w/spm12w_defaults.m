% spm12w_defaults
%
% spm12w defaults contains default parameters for preprocessing and glm 
% estimation. 
%
% To override these defaults the user has two options: 
%
% Option 1: If there are a set of parameters that the user wishes to override 
% for *ALL* data analyses simply copy them into the spm12w_userdefaults.m file 
% and place it somewhere in the the matlab path (ideally not in the global
% spm12w directory as that will override those values for all users sharing 
% that path). Any parameter in that file will override the spm12w defaults.
%
% Option 2: If there are parameters that the user wishes to override on a 
% study by study basis, place them at the bottom of your parameters file
% (for preprocessing) or glm file (for glm computation) and they will
% take precedence over spm12_defaults and spm12w_userdefaults.
%
% Option 3: You may also override the defaults using spm12w_userdefaults.m
% and override your own userdefaults on a study by study basis by placing 
% alternative default parameters at the bottom of your parameters file. 
%
% Note that when defaults are overriden, spm12w will indicate this in the
% command window output. If you DO NOT see these messages, then spm12w is
% not overridding the defaults. The most likely reason is that the
% spm12w_userdefaults.m file (option 1 above) is not in the matlab path or,
% if using options 2 or 3, that your parameters or glm file contains a typo. 
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: March, 2014 | Updated: June, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% spm12w naming conventions
def.fmri         = 'epi';            % Prefix for functional data
def.mri          = 'anat';           % Prefix for anatomical data
def.mristrip     = 'brain';          % Prefix for skull stirpped anatomical data
def.mnitok       = 'native';         % Current space of data
                                     % Fill in glm name at creation 
def.niceline     = [repmat('=',1,80)]; % A nice line for printing.

% spm12w logging defaults
% Logging level: 
%         1-All messages shown ([INFO],[WARNING],[EXCEPTION],[DEBUG])
%         2-Show only [INFO] & [WARNING] & [EXCEPTION]
%         3-Show only [INFO]
%         4-No messages (spm12w will run silently)
def.loglevel     = 1;
def.preplog      = sprintf('%s.log', p_.prep_name); % Default prepro logfile name
def.glmlog       = sprintf('%s.log', p_.glm_name);  % Default glm logfile name. 

% spm12w paths
% def.spm12wloc The dir where spm12w files are located 
% def.root:       The root dir for your study
% def.anadir:     The dir where user specific analyses are stored
% def.auxdir:     The dir where behavioral and other auxil data are stored.
% def.notedir:    The dir study related notes are stored
% def.prepdir:    The dir where the current preproc pipeline is stored
% def.qadir:      The dir where preprocessing qa output is stored
% def.rawdir:     The dir containing the "cleaned" raw nifti files
% def.scriptsdir: The dir where all user specific scripts are stored
% def.datadir:    The dir where preprocessing output is stored
% def.onsdir:     The dir where glm onsets are stored
% def.glmdir:     The dir where glm output will be saved
% def.rfxdir:     The dir where 2nd level rfx output will be saved
% def.roidir:     The dir where roi analyses output files will be saved
% def.imgfiledir: The dir where img files are location (e.g., bigmask, etc.)
% def.roimaskdir: The dir where image based masks are stored
% def.roispecdir: The dir where roi spec and variable files can be found
% def.archtok:    The dir name where prior analyses will be archived
def.spm12wloc  = fileparts(which('spm12w.m'));
def.root       = p_.study_dir;
def.anadir     = fullfile(def.root,'analysis',p_.username);
def.auxdir     = fullfile(def.root, 'auxil');
def.notesdir   = fullfile(def.root, 'notes');
def.prepdir    = fullfile(def.root, 'prep', p_.prep_name);
def.qadir      = fullfile(def.root, 'qa', p_.prep_name);
def.rawdir     = fullfile(def.root, 'raw', p_.sid);
def.scriptsdir = fullfile(def.root, 'scripts', p_.username);
def.datadir    = fullfile(def.prepdir, p_.sid);
def.onsdir     = fullfile(def.auxdir, 'onsets', p_.ons_dir);
def.glmdir     = fullfile(def.anadir, 'glm', p_.glm_name, p_.sid);
def.rfxdir     = fullfile(def.anadir, 'rfx', p_.rfx_name);
def.roidir     = fullfile(def.anadir, 'roi', p_.roi_name);
def.imgfiledir = fullfile(def.spm12wloc,'img_files');
def.roimaskdir = fullfile(def.spm12wloc,'roi_masks');
def.roispecdir = fullfile(def.auxdir,'roicsv');
def.archtok = 'archive'; 

% Utility function defaults
def.resample = 0; %Resampling for voxel extraction in spm12w_readnii
                   
% spm12w defaults for preprocessing
% Preprocessing defaults
def.trimvols = 0;           % Trim vols off front & back of session (default=0)
                            % can be vector [frontVols, backVols] i.e., [5,15].
                            % or matrix (sess X vols) [5,15; 2,12].                         
def.shuffle    = 1;         % Shuffle check
def.despike    = 0;         % Despiking
def.slicetime  = 1;         % Slicetime correction
def.realign    = 1;         % Motion correction
def.snr        = 1;         % Calculate snr
def.slices     = 1;         % Calculate slice noise

% Despiking defaults
def.cut        = [2.5,4];   % Afni defaults for despiking
def.c_order    = 'nvols/30';% Afni defaults for despiking
def.dpolort    = 2;         % Number of polynomials for despiking
def.least_sq   = 0;         % If 0 then do L1 regression for despiking
def.dilate     = 4;         % Dilates spm's "automask" by factor during despiking                          
                       
% Slicetime correction defaults
def.sformula   = 'philips'; % Slice order formula ('philips' or 'regular')
def.refslice   = 1;         % Reference slice for slice timing

% Realign defaults
def.realign_rtm = 1;     % Realign to 1st only (0) or 2pass(1st then sess mean) (1)
def.realign_quality = 1; % Quality versus speed trade-off. SPM12 default 0.9

% Normalization defaults
def.coreg2epi = 1;       % Coregister anatomy to mean epi image
def.ebiasreg   = 0.0001; % Intensity bias regularisation
def.ebiasfwhm  = 60;     % bias fwhm
def.etpm       = fullfile(spm('dir'),'tpm/TPM.nii'); % tissue prob maps
def.esamp      = 3;       % Sampling distance tradoff (smaller uses more data)


%spm8w leftovers, might need them
%def.mrf           = 2;        % Segmentation cleanup via MRF (0 | 1 | *2 | 3 | 4)
%def.warpreg       = 4;        % Warp regularisation (spm default:4|can use 0.4 for T1)                        
%def.sampling      = 3;        % Sampling distance (1|2|*3) < 3mm will slow computation 
%def.tissues       = 2;        % Tissue classes for template (1:GM|*2:GM+WM|3:GM+WM+CSF)
%def.dartelsmooth  = 6;        % Smooting kernel for DARTEL to MNI (GM,WM,EPI only)

% Normalization and write parameters
% Interpolation set to 7 bsplines gives a less blury normalized EPI.
% After smoothing that extra spatial resolution is lost.
% However, for studies with low smoothing, it might be worth bumping up
% the interp_r and interp_w and interp_c.
% r = realign, w = unwarp, c = coregister
def.wrap_r     = [0 0 0];     %SPM12/8 default [0 0 0] | spm2w was [0 1 0]
def.wrap_w     = [0 0 0];     %SPM12/8 default [0 0 0] | spm2w was [0 1 0]
def.interp_r   = 7;           %SPM12/8 default 4       | spm2w was 7
def.interp_w   = 7;           %SPM12/8 default 4       | spm2w was 7
def.interp_c   = 7;           %SPM12/8 default 4       | spm2w was 4
def.evoxsize   = [3 3 3];     %SPM12/8 default [2x2x2] | Don't resample
def.avoxsize   = [1 1 1];     %SPM12/8 default [2x2x2] | Don't resample
def.boundbox   = [-78 -114 -72;...
                   78 78 84]; %SPM12 default is [-78 -112 -70; 78 76 85]
                              %SPM8 default is [-78 -112 -50; 78 76 85]
                              %Ours accomodates with 3x3x3 and 2x2x2. 

% Slices noise defaults
def.noiseth = [5,15]; % Defaults for displaying slice noise figures. 
                      % May need tweaking for different scanners.

% Cleanup defaults
def.cleanup    = 1;  % 0 = keep all files (useful for debugging)
                     % 1 = delete all but the last stage of preprocessing (usually 's')
                     % 2 = delete all but the last 2 stages of preprocessing (usually 'sw')
                     % 3 = delete all but the last 2 stages of prerocessing and the raw epi
def.matkill    = 1;  % Delete intermediate mat files during cleanup
def.cleanupzip = 1;  % gzip final datasets (0 = no). 

% Brain mask
% Can use either bigmask in 1x1x1 or 3x3x3
def.mask = fullfile(def.imgfiledir,'bigmask_3x3x3.nii');

% GLM Model Specifications
def.runsplit    = 0;  % Seperate GLM per run?
def.design_only = 0;  % Design only (i.e., no data)
def.use_smooth  = 1;  % Run GLM on: 1 = smoothed, 2 = unsmoothed
def.duration    = 0;  % Event/Block Duration (same units as glm.time). 
                      % Dur files will override.
def.block_conv  = 1;  % Convolve blocks (0 for state/item) - 1=yes 0=no
def.polort      = 1;  % Order of polynomials: 0:Const|1:Linear|2:Quad|3:Cubic
def.outliers    = 1;  % Include outlier scans as nuissance? (use spm12w_art.m)
def.move        = 1;  % Include motion regressors?
def.nuissance   = 1;  % Include additional nuissance (run constants, constant)

% GLM Model Defaults
def.ons_ext    = 'txt';   % If no extension leave blank (ie: '';)   
def.durtime    = 'scans'; % Event duration in 'scans' (default) or 'sec' or 'ms'?
def.time       = 'scans'; % Onsets specified in 'scans' or 'secs' 
def.hpf        = Inf;     % HPF inf=no cutoff|otherwise cutoff in secs i.e. 128)                               
def.autocorr   = 'none';  % Autocorrelation correction (none | 'AR(1) + w')                                   
def.demean     = 0;       % Demean condition regressors as in SPM99?
def.orth       = 1;       % Disable w/i trial orth (0) or enable it (1, SPM default) 

% GLM HRF Specifications:
% OPTIONS:'hrf'
%         'hrf (with time derivative)'
%         'hrf (with time and dispersion derivatives)'
%         'Fourier set'
%         'Fourier set (Hanning)'
%         'Gamma functions'
%         'Finite Impulse Response'
def.hrf       = 'hrf';  %Set to Finite Impulse Response for MIXED designs.
def.hrfwindow = 20;     %For FIR: Length of HRF window 
def.hrfbasis  = 8;      %For FIR: Number of bins per window

% RFX Specification
% rfx_type : 'one-sample'
%            'two-sample' %not yet implemented
%            'anova1'     %not yet implemented
def.rfx_type  = 'one-sample';
def.rfx_im    = 1; % Implicit masking for NaN & 0 values at 2nd level.

% ROI Specification
def.roi_size  = 6;  % Default roi size is 6mm
def.spec_file = ''; % name of a csv file specifying rois
def.var_file  = ''; % name of a csv file with subject variables
