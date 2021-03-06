function spm12w_glm_compute(varargin)
% spm12w_glm_design('sub_id','glm_file')
%
% Inputs
% ------
% sid:      Subject ID of subject for glm computation (e.g., 's01')
%
% glm_file: File specifying the parameters for glm design and contrasts
%           (e.g., 'glm_tutorial.m'). If the full path is left unspecified,
%           spm12w_glm will look in the scripts directory. <optional>
%
% spm12w_glm_compute will gather onset files and specify a design matrix 
% which will be used to compute parameter estimates for a variety of designs.
% spm12w_glm_compute supports a variety of designs: event-related, block,
% regressor-only (e.g., PPI), state-item, and combinations of event, block
% and regressor designs. In addition, spm12w_glm_compute will generate
% outlier, movement and trends regressors to be included in the design
% matrix. 
%
% spm12w_glm_compute may also be used to specify a design without data if 
% the design_only parameter is set to 1 in the glm parameter file. In 
% this case the user must supply the tr, nses and nvols in the glm
% parameter file. Otherwise these parameters will be pulled from saved
% parameters in the prep directory specified by the user. 
%
% The first argument is a sid. The second argument is the name of a glm 
% parameter files (e.g., glm_tutorial.m) and is optional. If the parameter 
% file is unspecified, matlab will prompt the user to select the file.
%
% A pdf file (glmname.pdf) will be written to the analysis/username/glm_name 
% directory and contains a figure of the design matrix. 
%
% NB: Design matrices generated by spm12w_glm_compute build were tested against
% the same parameters in the SPM 12 gui and produced identical designs.
%
% Examples:
%
%       >>spm12w_glm_compute('sid','s01')
%       >>spm12w_glm_compute('sid','s01','glm_file', ...
%                           './scripts/username/glm_tutorial.m')
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: March, 2006 | Updated: June, 2017
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('sid','', 'glm_file','');
args = spm12w_args('nargs',2, 'defaults', args_defaults, 'arguments', varargin);

% Load glm parameters
glm = spm12w_getp('type','glm', 'sid',args.sid, 'para_file',args.glm_file);

% Setup directories for GLM analysis. 
spm12w_dirsetup('dirtype','glm','params',glm);

% Setup logfile
spm12w_logger('msg','setup_glm', 'level',glm.loglevel, 'params',glm)

% Goto glm directory
cd(glm.glmdir)

% If glm file has missing parameters, get them from the prep file. 
if isfield(glm, 'nses') && isfield(glm, 'nvols') && isfield(glm, 'tr') && isfield(glm,'nslice') && isfield(glm, 'refslice') && isfield(glm, 'sliceorder')
% Do nothing, estimation is go. 
elseif exist(fullfile(glm.datadir,[glm.prep_name,'.mat']),'file')
    load(fullfile(glm.datadir,[glm.prep_name,'.mat']))
    spm12w_logger('msg',sprintf(['Loading additional parameters from: ', ...
              '%s.mat'],glm.prep_name),'level',glm.loglevel)
    for pfield = {'nses','nvols','tr','nslice','refslice','sliceorder','ra','fmri','cleanupzip'};
        glm.(pfield{1}) = p.(pfield{1});
    end
    clear p % remove p structure for safety
else
    spm12w_logger('msg', ['[EXCEPTION] Required parameters (nses, ', ...
                  'nvols, tr, nsclice, refslice or sliceorder) are unspecified and a prep parameter file ', ...
                  'was not found.'], 'level',glm.loglevel)
    diary off
    error('Missing required parameters for model estimation (nses, nvols, tr, nslice, refslice or sliceorder)')     
end

% Show user the parameters we harvested.
for ses = 1:glm.nses
    spm12w_logger('msg', sprintf('Run:%d, nvols:%d, TR:%.2f', ses, ...
                  glm.nvols(ses), glm.tr(ses)), 'level', glm.loglevel);
end
spm12w_logger('msg', sprintf('nslice:%d, refslice:%d, sliceorder:%s', ...
              glm.nslice, glm.refslice, mat2str(glm.sliceorder)), ...
              'level', glm.loglevel);

% Adjust nses & nvols to match modeled runs.
if strcmp(glm.include_run,'all')
    glm.include_run = 1:glm.nses;
    idx_bool = ones(sum(glm.nvols),1); % keep all volumes
elseif max(glm.include_run > glm.nses)
    spm12w_logger('msg', sprintf(['[EXCEPTION] Included run (%d) exceeds ', ...
              'available runs (%d)'], max(glm.include_run),glm.nses), ...
              'level', glm.loglevel)
    error('Included run exceeds available runs...')
else
    % Generate an index of volumes to be modeled (for use in identifying onsets)
    idx_bool = zeros(sum(glm.nvols),1);
    vol = 1;
    for i_ses = 1:glm.nses
        for i_vols = 1:glm.nvols(i_ses)
            if any(i_ses==glm.include_run)
                idx_bool(vol) = 1;
            end
                vol = vol + 1;
        end
    end
    % Adjust nses & nvols to match modeled runs.
    glm.nses = length(glm.include_run);
    glm.nvols = glm.nvols(glm.include_run);
    % Assign TR to runs to be modeled
    glm.tr = glm.tr(glm.include_run);
end

% Check that the runs to be modeled all have same TR
if length(unique(glm.tr)) > 1
    spm12w_logger('msg', sprintf(['[EXCEPTION] Runs to be modeled have ', ...
              'different TRs (%s)'], sprintf('Run %.2f ',glm.tr)), ...
              'level', glm.loglevel)
    error('Modeled runs do not all have the same TR...')
else
    glm.tr = unique(glm.tr);
end

% Tell the user what we're up to...
msg = sprintf('GLM will be calculated on runs: %s', ...
               sprintf('%d(nvols=%d) ',[glm.include_run;glm.nvols]));
spm12w_logger('msg', msg, 'level', glm.loglevel)

% Build simple model from onsets and trends regressors.
% Note that constants are now specified outside of spm12w_glm_build
% because SPM demeans all user regressors and we want constants to stay
% as 1s and 0s so we specify them after the SPM design gets filled in by
% spm_fMRI_design.m. This is all because we're concatenating runs and is 
% equivalent to the tricks in spm12's spm_fmri_concatenate.m which also
% hacks on run constants after the design is filled in.
for mfield = {'events','blocks','regressors'};
    if ~isempty(glm.(mfield{1}))
        glm = spm12w_glm_build('type',mfield{1},'params',glm);
    end
end
for mfield = {'outliers','trends','move'}; %removed constants from build
    if glm.(mfield{1}) == 1
        glm = spm12w_glm_build('type',mfield{1},'params',glm);
    end
end   

% Adjust onsets, durations, parametrics and regressors for the included runs.
% If user is modeling fewer than the total number of sessions then 
% onsets (which always describe all runs) need to be adjusted by removing 
% the excluded run(s) and by altering the value of the onsets following the
% excluded run(s). 
if any(idx_bool == 0) && glm.runadjust == 1 % only adjust if necessary.
    spm12w_logger('msg',sprintf(['Adjusting onsets, durations,  ',...
                  'parametrics and/or regressors for the included runs: %s'], ...
                  mat2str(glm.include_run)), 'level',glm.loglevel)  
    for fname = fieldnames(glm.X_onsets)'
        if ismember(fname, glm.regressors)
            % If the onsets are regressors then different trick
            glm.X_onsets.(fname{1}).ons(idx_bool==0) = [];
        else
            % Onset type is block or events.
            % Create vectors describing every volume in the design
            ons_idx = zeros(length(idx_bool),1);
            % Mark 1 where onsets should be in the design     
            ons_idx(glm.X_onsets.(fname{1}).ons+1) = 1; % +1 for 0 to 1 indexing
            % Make new vectors for onsets and durations in same place
            ons = ons_idx;
            dur = ons_idx;
            % Put durations at same place as onsets
            dur(dur==1) = glm.X_onsets.(fname{1}).dur;
            % Remove all volumes for excluded runs based on idx_bool made above
            dur(idx_bool==0) = [];
            ons(idx_bool==0) = [];
            % Pop out all zeros from the new duration vector.
            dur(ons==0) = [];
            % Assign onsets and durations back to structure
            glm.X_onsets.(fname{1}).ons = find(ons==1)-1; % -1 for 1 to 0 indexing
            glm.X_onsets.(fname{1}).dur = dur;
            % Do same for parametrics
            if isfield(glm.X_onsets.(fname{1}),'P')
                for p_i = 1:length(glm.X_onsets.(fname{1}).P)
                    para = ons_idx;
                    para(para==1) = glm.X_onsets.(fname{1}).P(p_i).P;
                    para(idx_bool==0) = [];
                    para(ons==0) = [];
                    glm.X_onsets.(fname{1}).P(p_i).P = para;
                end            
            end
        end
    end
end

% Fill in SPM structure prior to design specification
glm.SPM = spm12w_getspmstruct('type','glm','params',glm);

% Do seperate steps if design_only  == 1 or 0
if glm.design_only == 0
    spm12w_logger('msg','Generating SPM design matrix...','level',glm.loglevel)
    % Print files being used
    for file_i = 1:size(glm.SPM.xY.Pfiles,1)
        epifile = deblank(glm.SPM.xY.Pfiles(file_i,:));
        spm12w_logger('msg',sprintf('GLM will use file: %s', epifile), ...
                      'level',glm.loglevel)      
    end
    
    % Generate design using SPM structure (no figure for now)
    glm.SPM = spm_fmri_spm_ui(glm.SPM);
    
    % Append constants.
    runcons = [];
    con_names = {};
    for ses_i = 1:glm.nses
        spm12w_logger('msg', sprintf(['[DEBUG] Generating run constants ', ...
                      'for run: %d'], glm.include_run(ses_i)),'level',glm.loglevel)
        runcons = blkdiag(runcons,ones(glm.nvols(ses_i),1));
        con_names{ses_i} = sprintf('Sn(%i) constant',ses_i);        
    end
    
    % Place constants in design and adjust the iB field.
    glm.SPM.xX.X    = [glm.SPM.xX.X(:,1:end-1) runcons];
    glm.SPM.xX.iB   = glm.SPM.xX.iB:(glm.SPM.xX.iB+size(runcons,2)-1);
    glm.SPM.xX.name = {glm.SPM.xX.name{1:end-1} con_names{:}}; 
    
    % Now adjust the iG (nuissance) field for the user entered regressors.
    % We have to do this manually as SPM as no way of knowing if user specified
    % regressors (i.e., linear trends and motion) are nuissance or not.
    % Find location of trends, move and outliers
    nuiss_idx = [];
    for nuissance = {'outliers','trends','move'}
        nuiss_str = sprintf('Sn(1) r-%s',nuissance{1});
        nuiss_idx = [nuiss_idx, find(strcmp(glm.SPM.xX.name, nuiss_str))];
    end
    % Let user know what we found
    if ~isempty(nuiss_idx)
        spm12w_logger('msg', sprintf(['[DEBUG] Setting the following design ', ...
                      'columns to nuissance: %s'], mat2str(nuiss_idx)), ...
                      'level',glm.loglevel)  
        glm.SPM.xX.iG = nuiss_idx;
        glm.SPM.xX.iC(nuiss_idx) = [];
    end
              
    % Demean design if user requested
    if glm.demean == 1
        reg_interest = size(glm.SPM.Sess.U,2);
        for i = 1:size(glm.SPM.Sess.U,2)
            for ii = 1:size(glm.SPM.Sess.U(i).P,2)
                if strfind(glm.SPM.Sess.U(i).P(ii).name,'other')
                    reg_interest = reg_interest + 1;
                end
            end
        end
        for i = 1:reg_interest
            glm.SPM.xX.X(:,i)=spm_detrend(glm.SPM.xX.X(:,i),0);
        end
    end
       
    % Modify SPM structure for runwise HRF (from spm_fmri_concatenate.m)
    if glm.hpf ~= Inf
        spm12w_logger('msg',sprintf(['Including high pass filtering ',...
                  'in design using hpf: %d'], glm.hpf), 'level',glm.loglevel)  
        s = cumsum([0 glm.nvols]); 
        for i=1:numel(glm.nvols)
            K(i) = struct('HParam', glm.SPM.xX.K(1).HParam,...
                          'row',    s(i) + (1:glm.nvols(i)),...
                          'RT',     glm.SPM.xY.RT);
        end
        glm.SPM.xX.K = spm_filter(K);
    end
    
    % Apply runwsie autocorrelation correction
    switch lower(glm.SPM.xVi.form)
        case {'ar(1)','ar(0.2)'}
            spm12w_logger('msg',sprintf(['Serial autocorrelation ',...
                  'correction type: %s'], glm.autocorr), 'level',glm.loglevel)  
            glm.SPM.xVi.Vi   = spm_Ce(glm.nvols,0.2);
            glm.SPM.xVi.form = 'AR(0.2)';
        case {'i.i.d', 'none'}
        otherwise
            warning('Unhandled temporal non-sphericity.');
    end
      
    % Specify mask
    [~,maskname] = fileparts(glm.mask);
    glm.SPM.xM.VM = spm_vol(glm.mask);
    glm.SPM.xM.T  = [];
    glm.SPM.xM.TH = ones(size(glm.SPM.xM.TH))*(-Inf);
    glm.SPM.xM.I  = 0;
    glm.SPM.xM.xs = struct('Masking', sprintf('explicit masking only - using %s',maskname));
    spm12w_logger('msg',sprintf(['The SPM.mat file has been modified to ',...
                  'use mask: %s'], maskname), 'level',glm.loglevel)                               
              
    % Estimate model
    spm12w_logger('msg',sprintf('Estimating parameters for model: %s', ...
                  glm.glm_name), 'level',glm.loglevel) 
    glm.SPM = spm_spm(glm.SPM);
    
    % Generate an empty figure to catch the estimated SPM design
    F = spm_figure('CreateWin','Graphics', 'spm12w glm', 'off'); 
    fname = reshape(cellstr(glm.SPM.xY.P),size(glm.SPM.xY.VY));
    spm_DesRep('DesMtx',glm.SPM.xX,fname,glm.SPM.xsDes)
    
    % Hide F right away (unlike preprocessing, spm_fmri_spm_ui will
    % set a hidden figure to visible). 
    set(F,'visible','off');
    % Print looks better using opengl for design matrices. Keep an eye on
    % this in case it fails on other platforms. opengl throws running in parfor
    % but completed anyway.
    print(F, 'glm.ps', '-dpsc2','-opengl','-append','-noui')   
       
    % Add effects of interest contrast (nice to have and necessary for
    % adjusting timecourses during VOI extraction for PPI). 
    Fcname = 'effects of interest';
    iX0      = [glm.SPM.xX.iG glm.SPM.xX.iB ]; % span nuissance + constants
    glm.SPM.xCon = spm_FcUtil('Set',Fcname,'F','iX0',iX0,glm.SPM.xX.xKXs);
    % Save to SPM structure
    SPM = glm.SPM;
    save('SPM.mat','SPM')       
    % Write model residuals adjusting for the effects of interest contrast
    if glm.writeresid
        % Todo allow flexbile contrast selection for
        % adjustment. For now I'm hardcoding to the effects of interest 
        % contrast which will remove everything but the task effects.
        % (i.e., demeans run constants, motion, linear trends, etc.)
        spm_write_residuals(glm.SPM,1);
        % Create list of files to concat (using sprintf trick)
        resfiles = strsplit(sprintf('Res_%04d.nii ',1:glm.SPM.nscan))';
        % spm_file_merge expects a char array so convert the cell array
        resfiles = char(resfiles(1:glm.SPM.nscan));     
        spm_file_merge(resfiles,glm.residname,0);
        gzip(glm.residname)
        % Remove the unzipped resdiuals file
        spm_unlink(glm.residname)
        % Remove the temp mat file that spm created 
        spm_unlink(strrep(glm.residname,'.nii','.mat'))
        % Remove the 3d residual files
        for delfile = cellstr(resfiles)'        
            spm_unlink(delfile{1})
        end
    end
    % Remove the unzipped files if they were previously zipped
    if glm.cleanupzip == 1 
        for file_i = 1:size(glm.SPM.xY.Pfiles,1)
            epifile = deblank(glm.SPM.xY.Pfiles(file_i,:));
            delete(epifile)
        end    
    end    
    % Convert multipage ps file to pdf
    spm12w_ps2pdf('ps_file',fullfile(glm.glmdir,'glm.ps'),...
                  'pdf_file',fullfile(glm.glmdir,[glm.glm_name,'.pdf']))
    
    % Save glm parameters structure              
    save(sprintf('%s.mat',glm.glm_name),'glm')    
    
    % Set final words
    % Determine design type
    if ~isempty(glm.events) && isempty(glm.blocks) && isempty(glm.regressors)
        destype = 'event-related';
    elseif isempty(glm.events) && ~isempty(glm.blocks) && isempty(glm.regressors)
        destype = 'block';
    elseif isempty(glm.events) && isempty(glm.blocks) && ~isempty(glm.regressors)
        destype = 'regressor only';
    elseif ~isempty(glm.events) || ~isempty(glm.blocks) && isempty(glm.regressors)
        destype = 'Mixed Design (e.g., State-Item)';  
    elseif ~isempty(glm.events) || ~isempty(glm.blocks) && ~isempty(glm.regressors)
        destype = 'unknown hybrid regressor or PPI';  
    end       
    if ~isempty(glm.parametrics)
        destype = sprintf('%s | parametric modulation', destype);
    end
        
    msglist{1} = glm.niceline;
    msglist{2} = sprintf('GLM specification complete on subject: %s',glm.sid);
    msglist{3} = sprintf('Design Type: %s',destype);
    msglist{4} = sprintf('Summary    : TR:%.2f, Runs:%d, Nvols:%s, Nslices:%d',unique(glm.tr), glm.nses,mat2str(glm.nvols), glm.nslice);
    msglist{5} = sprintf('Summary    : hfr:%s, hpf:%d, autocorr:%s, time:%s',glm.hrf, glm.hpf,glm.autocorr, glm.time);
    msglist{6} = sprintf('Summary    : microtime resolution:%d, microtime onset :%d',glm.tbins, glm.tref);
    msglist{7} = sprintf('Summary    : Regressors: effects(%d), nuissance(%d), constants(%d)',numel(glm.SPM.xX.iC), numel(glm.SPM.xX.iG), numel(glm.SPM.xX.iB));
    msglist{8} = sprintf('Log file   : %s', fullfile(glm.glmdir, glm.glmlog));
    msglist{9} = sprintf('Figures    : %s', fullfile(glm.glmdir,[glm.glm_name,'.pdf']));
    msglist{10} = sprintf('Parameters : %s', fullfile(glm.glmdir,[glm.glm_name,'.mat']));
else   
    spm12w_logger('msg','Generating SPM design matrix (design only)...','level',glm.loglevel)
    
    % Generate design using SPM structure.
    glm.SPM = spm_fmri_spm_ui(glm.SPM);
    
    % Append constants.
    runcons = [];
    con_names = {};
    for ses_i = 1:glm.nses
        spm12w_logger('msg', sprintf(['[DEBUG] Generating run constants ', ...
                      'for run: %d'], glm.include_run(ses_i)),'level',glm.loglevel)
        runcons = blkdiag(runcons,ones(glm.nvols(ses_i),1));
        con_names{ses_i} = sprintf('Sn(%i) constant',ses_i);        
    end
    
    % Place constants in design and adjust the iB field.
    glm.SPM.xX.X    = [glm.SPM.xX.X(:,1:end-1) runcons];
    glm.SPM.xX.iB   = glm.SPM.xX.iB:(glm.SPM.xX.iB+size(runcons,2)-1);
    glm.SPM.xX.name = {glm.SPM.xX.name{1:end-1} con_names{:}}; 
    
    % Now adjust the iG (nuissance) field for the user entered regressors.
    % We have to do this manually as SPM as no way of knowing if user specified
    % regressors (i.e., linear trends and motion) are nuissance or not.
    % Find location of trends, move and outliers
    nuiss_idx = [];
    for nuissance = {'outliers','trends','move'}
        nuiss_str = sprintf('Sn(1) r-%s',nuissance{1});
        nuiss_idx = [nuiss_idx, find(strcmp(glm.SPM.xX.name, nuiss_str))];
    end    
    % Let user know what we found
    if ~isempty(nuiss_idx)
        spm12w_logger('msg', sprintf(['[DEBUG] Setting the following design ', ...
                      'columns to nuissance: %s'], mat2str(nuiss_idx)), ...
                      'level',glm.loglevel)  
        glm.SPM.xX.iG = nuiss_idx;
        glm.SPM.xX.iC(nuiss_idx) = [];
    end
       
    % Generate an empty figure to catch the SPM design
    F = spm_figure('CreateWin','Graphics', 'spm12w glm', 'off'); 
    fname = ''; %no files if design only
    spm_DesRep('DesMtx',glm.SPM.xX,fname,glm.SPM.xsDes)
    
    % Hide F right away (unlike preprocessing, spm_fmri_spm_ui will
    % set a hidden figure to visible). 
    set(F,'visible','off');
    % Print looks better using opengl for design matrices. Keep an eye on
    % this in case it fails on other platforms. opengl throws running in parfor
    % but completed anyway.
    print(F, 'glm.ps', '-dpsc2','-opengl','-append','-noui')       
    
    % Save to SPM structure
    SPM = glm.SPM;
    save('SPM.mat','SPM')   
    
    % Save glm parameters structure              
    save(sprintf('%s.mat',glm.glm_name),'glm')  
    
    % Set final words (no figure, design only)
    msglist{1} = glm.niceline;
    msglist{2} = sprintf('GLM specification complete on subject: %s',glm.sid);
    msglist{3} = sprintf('Log file   : %s', fullfile(glm.glmdir, glm.glmlog));
    msglist{4} = sprintf('Parameters : %s', fullfile(glm.glmdir,[glm.glm_name,'.mat']));
end

% Print final words
for msg = msglist
    spm12w_logger('msg',msg{1},'level',glm.loglevel)
end

% Close hidden figure (try because in some cases it might already be closed)
try
    F = spm_figure('FindWin','Graphics');
    close(F)
end

% Close log and return to studydir.
diary off; 
cd(glm.study_dir)