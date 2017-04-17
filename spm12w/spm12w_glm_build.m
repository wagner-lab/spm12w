function glm = spm12w_glm_build(varargin)
% spm12w_glm_build('type','','params','')
%
% Inputs
% ------
% type   : Type of glm model to build. Options are 'events', 'blocks', 
%          'regressors', 'outliers', 'move', 'trends' (default = 'events').
%
% params : A structure of glm parameters (usually loaded from a glm
%          parameters file via spm12w_getp. 
%
% spm12w_glm_build will gather onset files and output a design matrix. In
% practice spm12w_glm_build is called by spm12w_glm_compute multiple times
% to build a complete design matrix (i.e., once for events, once for
% including outliers and again for including movement parameters). 
%
% The first argument is a model type. The second argument is a glm
% parameter structure. Model types are:
%
% The following model types are task conditions and may be confovled with HRF
% 'events'     : Standard rapid event related design (may have durations)
% 'blocks'     : Standard block design with option to convolve or not with HRF
% 'regressors' : Regressor designs (e.g., PPI) usually not convolved with HRF
%
% The following model types are trends variables and never convovled with HRF
% 'outliers'  : Outlier volumes (for de-weighting bad scans)
% 'move'      : Movement/realignment parameters (optionally add 1st derivative)
% 'trends'    : Standard trends regressors
% 'constants' : Standard run constant regressors for n-1 runs. {DEPRECIATED}
%
% Examples:
%
%       >>spm12w_glm_build('type','events','params',glm)
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: March, 2006 | Updated: March, 2017
% # TODO: Figure out how to disable spm orthogonalization of regressors as
% #       spm12 has added the ability to do so. 
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('type','events', 'params','');
args = spm12w_args('nargs',4, 'defaults', args_defaults, 'arguments', varargin);

% Assign params to glm for cleaner code
glm = args.params;

% Check if parametrics exist, as users may often not include the field in their 
% params file, we'll create it if it doesn't exist so as not to crash
% the model generation below
if ~isfield(glm,'parametrics')
    glm.parametrics = {};
end

% Load onset files and create a matrix of model predictors
% Column 1 = Onsets, Column 2 = Durations, Column 3...n = parametrics.
if ismember(args.type, {'events','blocks','regressors'})
    for onsname = glm.(args.type) 
        % Check for sid specific onsets otherwise global
        sidfname = sprintf('%s_%s.%s', glm.sid, onsname{1}, glm.ons_ext);
        onsfname = sprintf('%s.%s', onsname{1}, glm.ons_ext);
        if exist(fullfile(glm.onsdir, sidfname),'file')
            spm12w_logger('msg', sprintf('[DEBUG] Loading %s onsets for %s (file:%s)', ...
                          args.type, onsname{1},sidfname),'level',glm.loglevel)
            onsets.(onsname{1}).ons = spm_load(fullfile(glm.onsdir, sidfname));
        elseif exist(fullfile(glm.onsdir,onsfname),'file') 
            spm12w_logger('msg', sprintf('[DEBUG] Loading %s onsets for %s (file:%s)', ...
                          args.type, onsname{1},onsfname),'level',glm.loglevel)
            onsets.(onsname{1}).ons = spm_load(fullfile(glm.onsdir, onsfname));
        else
            spm12w_logger('msg', sprintf(['[EXCEPTION] Cannot find the ', ...
                  'onset file: %s for condition: %s... are you sure it ', ...
                  'exists?'], onsfname, onsname{1}), 'level',glm.loglevel)
            diary off
            error('Cannot find the onset file: %s for condition: %s...',...
                  onsfname, onsname{1})   
        end
        % Check for sid specific durations otherwise global
        sidfname = sprintf('%s_%s_dur.%s', glm.sid, onsname{1}, glm.ons_ext);
        onsfname = sprintf('%s_dur.%s', onsname{1}, glm.ons_ext);
        if exist(fullfile(glm.onsdir, sidfname),'file')
            spm12w_logger('msg', sprintf('[DEBUG] Loading %s durations for %s (file:%s)', ...
                          args.type, onsname{1},sidfname),'level',glm.loglevel)
            onsets.(onsname{1}).dur = spm_load(fullfile(glm.onsdir, sidfname));
        elseif exist(fullfile(glm.onsdir,onsfname),'file') 
            spm12w_logger('msg', sprintf('[DEBUG] Loading %s durations for %s (file:%s)', ...
                          args.type, onsname{1},onsfname),'level',glm.loglevel)
            onsets.(onsname{1}).dur = spm_load(fullfile(glm.onsdir, onsfname));
        else
            onsets.(onsname{1}).dur = repmat(glm.duration, ...
                                              length(onsets.(onsname{1}).ons),1);
        end
        % Adjust durations in case of time/durtime mismatch.
        if strcmp(glm.time,'scans') && strcmp(glm.durtime, 'sec')
            spm12w_logger('msg', sprintf(['[WARNING] Onset time is scans ',...
                          'but duration time is sec. Dividing durations by ',...
                          '%1.2f'], glm.tr), 'level',glm.loglevel)
            onsets.(onsname{1}).dur = onsets.(onsname{1}).dur./glm.tr; 
        elseif strcmp(glm.time,'scans') && strcmp(glm.durtime, 'ms')
            spm12w_logger('msg', sprintf(['[WARNING] Onset time is scans ',...
                          'but duration time is ms. Dividing durations by ',...
                          '%1.2f'], glm.tr*1000), 'level',glm.loglevel)
            onsets.(onsname{1}).dur = onsets.(onsname{1}).dur./(glm.tr*1000);        
        elseif strcmp(glm.time,'sec') && ~strcmp(glm.durtime, 'sec')
            spm12w_logger('msg', sprintf(['[WARNING] Onset time is sec ',...
                          'but duration time is scans. Multiplying durations by ',...
                          '%1.2f'], glm.tr), 'level',glm.loglevel)
            onsets.(onsname{1}).dur = onsets.(onsname{1}).dur.*glm.tr; 
        end
        % Check for special keyword to grab all the parametrics
        if size(glm.parametrics,2) == 1 && strcmp(glm.parametrics{1},'allthethings')
            % Check for sid specific parametrics otherwise global
            % As we don't know exactly how many parametrics, we do a dir. 
            sidfname = sprintf('%s_%sx*.%s', glm.sid, onsname{1}, glm.ons_ext);
            onsfname = sprintf('%sx*.%s', onsname{1}, glm.ons_ext);
            parlist = dir(fullfile(glm.onsdir,sidfname));
            if isempty(parlist)
                % If parlist based on sid is empty then try looking for onsfname
                parlist = dir(fullfile(glm.onsdir,onsfname));
            end
        elseif ~isempty(glm.parametrics)
            % Not elegant, but iterate through all parametrics and see if
            % they apply to current onset condition
            parlist = '';
            for para = glm.parametrics
                if strcmp([onsname{1},'x'],para{1}(1:length(onsname{1})+1))
                    % Check for sid specific parametrics otherwise global
                    % As we don't know exactly how many parametrics, we do a dir. 
                    sidfname = sprintf('%s_%s.%s', glm.sid, para{1}, glm.ons_ext);
                    onsfname = sprintf('%s.%s', para{1}, glm.ons_ext);
                    if exist(fullfile(glm.onsdir,sidfname),'file') == 2
                        parlist(end+1,1).name = sidfname;
                    elseif exist(fullfile(glm.onsdir,onsfname),'file') == 2
                        parlist(end+1,1).name = onsfname;
                    end   
                end
            end        
        else
            % No parametrics (empty list)
            parlist = '';      
        end
        if ~isempty(parlist)
            for par_i = 1:size(parlist,1)
                spm12w_logger('msg', sprintf('[DEBUG] Loading %s parametrics for %s (file:%s)', ...
                          args.type, onsname{1},parlist(par_i).name),'level',glm.loglevel)                
                parname = parlist(par_i).name;
                [~,parname,~] = fileparts(parname(strfind(parname,[onsname{1},'x'])+length(onsname{1})+1:end));
                onsets.(onsname{1}).P(par_i).name = parname;
                onsets.(onsname{1}).P(par_i).P = spm_load(fullfile(glm.onsdir, parlist(par_i).name)); 
            end            
        end
    end
    % Assign to glm structure
    for onsfield = fields(onsets)'
        glm.X_onsets.(onsfield{1}) = onsets.(onsfield{1});
    end
end

% Create outlier regressors from outlier files (one column per regressor)
if strcmp(args.type, 'outliers')
    outvols = [];
    for ses_i = 1:glm.nses
        outfile = sprintf('outliers_r%02d.%s',glm.include_run(ses_i),  glm.ons_ext);
        if exist(fullfile(glm.datadir,outfile),'file')
            spm12w_logger('msg', sprintf('[DEBUG] Loading outliers for run: %d (file:%s)', ...
                          glm.include_run(ses_i), outfile),'level',glm.loglevel)
            out_tmp = load(fullfile(glm.datadir,outfile));
            if ses_i > 1
                %Adjust to take into account session durations
                out_tmp = out_tmp + sum(glm.nvols(1:ses_i-1));
            end
            outvols = [outvols;out_tmp];
        end
    end
    if ~isempty(outvols) 
        % create a nvol x length of outliers matrix.
        outliers = zeros(sum(glm.nvols),length(outvols));
        for out_i = 1:length(outvols)
            outliers(outvols(out_i),out_i) = 1;
        end
        % Assign to glm structure
        glm.X_outliers = outliers;
    end
end

% Create movement regressors from ra variable in glm structure
% (we don't need to use txt file since we save the realignments to the p strct) 
if strcmp(args.type, 'move')
    move = [];
    for ses_i = 1:glm.nses
        spm12w_logger('msg', sprintf('[DEBUG] Loading realignment parameters for run: %d', ...
                      glm.include_run(ses_i)),'level',glm.loglevel)
        move{ses_i} = glm.ra{glm.include_run(ses_i)};
    end
    glm.X_move = cat(1,move{:});
end

% Create trends regressors
if strcmp(args.type, 'trends')
    % Make empty matrix to account for all the polys
    polys = zeros(sum(glm.nvols), glm.polort*glm.nses);
    row_n = 1;
    col_m = 1;
    % Iterate through sessions and columns assigning polys.
    for ses_i = 1:glm.nses
        spm12w_logger('msg', sprintf('[DEBUG] Generating polynomial trends (polort:%d) for run: %d', ...
                      glm.polort,glm.include_run(ses_i)),'level',glm.loglevel)
        for ii = 1:glm.polort
            polys(row_n:row_n+glm.nvols(ses_i)-1, col_m) = linspace(-1,1,glm.nvols(ses_i)).^ii';
            col_m = col_m + 1;
        end
        row_n = row_n + glm.nvols(ses_i);
    end
    glm.X_trends = polys;
end

% Create run constants
% This is depreciated March 2017 as we now add constants after SPM structure is
% filled in or else spm_fMRI_design.m will demean these. I'm leaving this code
% here in case we decide to revert.
% if strcmp(args.type, 'constants')
%     % Make empty matrix to account for all the polys
%     runcons = zeros(sum(glm.nvols), glm.nses-1);
%     row_n = 1;
%     col_m = 1;
%     for ses_i = 1:glm.nses
%         spm12w_logger('msg', sprintf('[DEBUG] Generating run constants (n-1 runs) for run: %d', ...
%                       glm.include_run(ses_i)),'level',glm.loglevel)
%         %runcons(row_n:row_n+glm.nvols(ses_i)-1, col_m) = linspace(-1,1,glm.nvols(ses_i)).^0';
%         runcons(row_n:row_n+glm.nvols(ses_i)-1, col_m) = ones(glm.nvols(ses_i),1);
%         col_m = col_m + 1;
%         row_n = row_n + glm.nvols(ses_i);
%     end
%     glm.X_constants = runcons;
% end