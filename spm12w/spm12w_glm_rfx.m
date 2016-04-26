function spm12w_glm_rfx(varargin)
% spm12w_glm_rfx('glm_file','sids','rfx_dir')
%
% Inputs
% ------
% glm_file: File specifying the parameters for glm design and contrasts
%           (e.g., 'glm_tutorial.m'). If the full path is left unspecified,
%           spm12w_glm will look in the scripts directory.
%
% sids:     A cell array of Subject IDs for random effects analysis. If
%           left unspecified, a dialog box will appear asking the user to
%           select subjects. Use the keyword 'allsids' and all subjects in
%           the specified glm directory will be used for rfx analysis.
%
% rfx_dir : Name of the directory to which results will be written. If left
%           unpecified, the directory name will be taken from the
%           glm.rfx_name variable in the glm parameters file. If that variable 
%           is absent the directory will have the same name as the GLM that was 
%           used for generating contrasts. <optional> 
%
% spm12w_glm_rfx will gather contrast files from a prior GLM for which contrasts
% have been generated and will run a random effects analysis on the resulting
% contrast files. At the moment, spm12w_glm_rfx supports only one-sample t-tests
% across subjects. In the future, this may expand to other potential
% models (e.g., anova, 2-sample t-test).
%
% The first argument is the name of a glm parameters file (e.g., glm_turoial.m).
% The second argument (optional) is a cell array of sids (can be left blank to 
% manually select them).
%
% The results of the random effects analysis will be output to a folder
% with name specified by rfx_name argument (optional).
%
% Examples:
%
%       >>spm12w_glm_rfx 
%       >>spm12w_glm_rfx('glm_file', './scripts/username/glm_tutorial.m', ...
%                        'sids', {'s01','s02','s03'})
%       >>spm12w_glm_rfx('glm_file', './scripts/username/glm_tutorial.m', ...
%                        'sids', {'allsids'}, ...
%                        'rfx_dir', 'rfx_tutorial')
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: March, 2006 | Updated: September, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('glm_file','', 'sids','','rfx_dir','');
args = spm12w_args('nargs',0, 'defaults', args_defaults, 'arguments', varargin);

% Load glm parameters
rfx = spm12w_getp('type','glm', 'para_file',args.glm_file);

% If rfx_dir argument was not provided, get rfx_dir from glm parameters
if ~isempty(args.rfx_dir) && isempty(strfind(args.rfx_dir,fullfile(rfx.anadir,'rfx')))
    rfx.rfxdir = fullfile(rfx.anadir, 'rfx', args.rfx_dir);
elseif ~isempty(args.rfx_dir)
    rfx.rfxdir = args.rfx_dir;  
end

% Check for cell in case user provided allsids as string.
if ~iscell(args.sids) && ~isempty(args.sids)
    args.sids = cellstr(args.sids);
end

% If sids argument was not provided, open dialog window to get sids.
% If sids argument contained the keyword 'allsids', then get all sids.
% Since we should only do rfx on computed glms, let's look in rfx.glmdir.
if isempty(args.sids)
    args.sids = spm12w_getsid(rfx.glmdir);
elseif numel(args.sids) == 1 && strcmp(args.sids,'allsids')
    sids = dir(fullfile(rfx.glmdir,'s*'));
    args.sids = {sids.name};
end

% Setup directories for RFX analysis. Archive prior spmT file and spm.mat.
spm12w_dirsetup('dirtype','rfx','params',rfx);

% Perform rfx analysis on the specified conditions
for rfxcondir = rfx.rfx_conds
    spm12w_logger('msg',sprintf(['Performing 2nd level rfx analysis (%s) ',...
                  'on contrast: %s'], rfx.rfx_type, rfxcondir{1}), ...
                  'level',rfx.loglevel)
    % Add the output dir for this rfx analysis to rfx structure and cd to it
    rfx.rfxcondir = fullfile(rfx.rfxdir,rfxcondir{1});
    cd(rfx.rfxcondir);
    % RFX analysis setup. 
    % Add file variables to rfx structure.
    % spm12w_getspmstruct will use this to generate the SPM structure. 
    rfx.rfxconfiles = {};
    for sid = args.sids
        % Load the SPM file for the GLM (do this subjectwise to be safe)
        SPM_ = load(fullfile(rfx.glmdir,sid{1},'SPM.mat'));
        % Find the index of the rfx contrast in the SPM.xCon
        conidx = find(strcmp({SPM_.SPM.xCon.name},rfxcondir{1}));
        % Use the index to get the filename of the con file for that sid
        rfx.rfxconfiles{end+1,1} = fullfile(rfx.glmdir,sid{1},...
                                          SPM_.SPM.xCon(conidx).Vcon.fname);
    end
    % Adjust paths to include spm12/config dir
    addconfig = fullfile(fileparts(which('spm.m')),'config');
    addpath(addconfig)
    % Generate the SPM structure for rfx
    rfx.SPM = spm12w_getspmstruct('type','rfx','params',rfx);  
    % Estimate the rfx model (using spm_run_fmri_est.m)
    job_est.spmmat = {fullfile(rfx.rfxcondir,'SPM.mat')};
    job_est.write_residuals = 0;
    job_est.method = struct('Classical',1);
    spm12w_logger('msg',sprintf(['Estimating 2nd level rfx model ',...
                  'on contrast: %s'], rfxcondir{1}), 'level',rfx.loglevel)
    spm_run_fmri_est(job_est);
    % Generate contrasts for the model
    job_con.spmmat = {fullfile(rfx.rfxcondir,'SPM.mat')};
    job_con.consess{1}.tcon = struct('name',rfxcondir{1},'weights',1,...
                                     'sessrep','none');
    job_con.delete = 0;
    spm12w_logger('msg',sprintf(['Generating statistics for 2nd level rfx ',...
              'model on contrast: %s'], rfxcondir{1}), 'level',rfx.loglevel)
    spm_run_con(job_con);    
    % rm spm12/config from path
    rmpath(addconfig)
    % Save the rfx structure
    save(sprintf('%s.mat',rfx.rfx_name),'rfx');
    % Remove the newly created fields in preperation for next iteration.
    rfx = rmfield(rfx,{'rfxcondir','rfxconfiles'});
end

% Set final words
msglist{1} = rfx.niceline;
msglist{2} = sprintf('2nd level rfx analysis (%s) complete...', rfx.rfx_name);
for rfxcondir = rfx.rfx_conds
    msglist{end+1} = sprintf('rfx analysis (%s) computed for contrast: %s', ...
                             rfx.rfx_type, rfxcondir{1});
    msglist{end+1} = sprintf('Parameters saved to : %s', ...
                     fullfile(rfx.rfx_name,rfxcondir{1},[rfx.rfx_name,'.mat']));
end

% Print final words
for msg = msglist
    spm12w_logger('msg',msg{1},'level',rfx.loglevel)
end

% return to studydir 
cd(rfx.study_dir)

