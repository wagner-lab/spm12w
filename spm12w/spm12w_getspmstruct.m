function SPM = spm12w_getspmstruct(varargin)
% spm12w_getspmstrcut(type, params)
%
% Input
% -----
% type   : Type of spm structure to build (default = 'glm'). 
%
% params : A structure of parameters (i.e., p, glm, etc.). 
%
% Takes a parameters structure and a spm structure type and fills in the spm
% structure based on the parameters structure. Type can be:
%
%   'glm' : spm12w_getspmstruct will create and fill in the spm structure 
%           appropriate for a 1st level GLM analysis. 
%
%   'rfx' : spm12w_getspmstruct will create and fill in the spm job structure
%           appropriate for subsequent generation of a full SPM structure. 
%           Unlike with 'glm', this script will not generate the full spm
%           structure, opting instead to let spm12w_glm_rfx.m generate the
%           appropriate SPM.mat file for 2nd level analysis (via calls to
%           spm_run_factorial.design.m). The reason for this is that, unlike
%           1st level analyses where we depart from the SPM norms, our 2nd 
%           level analyses are standard SPM and can rely on the SPM machienry
%           to a higher degree than 1st level analyses. Finally, the type of rfx
%           analysis is specified by the rfx_type parameter in the glm 
%           parameters file.
%
% Examples:
%
%   >> SPM = spm12w_getspmstruct('type, 'glm', 'params', glm)
%   >> JOB = spm12w_getspmstruct('type, 'rfx', 'params', rfx)
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: January 2015 | Updated: March, 2017
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('type','glm', 'params','');
args = spm12w_args('nargs',4, 'defaults', args_defaults, 'arguments', varargin);

% Assign params to p_spm (i.e, parameters spm) for cleaner code. 
p_spm = args.params;

% Setup SPM structure
switch args.type
    case 'glm'
        % Setup SPM Structure for 1st level GLM
        spm12w_logger('msg',sprintf(['[DEBUG] Setting up SPM structure for ',...
                      'type:%s and subject:%s'], args.type, p_spm.sid),...
                      'level',p_spm.loglevel)
        % User-defined parameters for this analysis
        SPM.nscan          = sum(p_spm.nvols);        % number of scans for each of nsess sessions
        SPM.xY.RT          = p_spm.tr;                % experiment TR in seconds
        SPM.xGX.iGXcalc    = 'None';                  % global normalization: OPTIONS:'Scaling'|'None'
        SPM.xX.K.HParam    = p_spm.hpf;               % high-pass filter cutoff (secs) [Inf = no filtering] 
        SPM.xVi.form       = p_spm.autocorr;          % intrinsic autocorrelations: OPTIONS: 'none'|'AR(1)' 
        % basis functions and timing parameters
        SPM.xBF.name       = p_spm.hrf;     % OPTIONS:'hrf'
                                            %         'hrf (with time derivative)'
                                            %         'hrf (with time and dispersion derivatives)'
                                            %         'Fourier set'
                                            %         'Fourier set (Hanning)'
                                            %         'Gamma functions'
                                            %         'Finite Impulse Response'

        SPM.xBF.T          = p_spm.tbins;   % number of time bins per scan
        SPM.xBF.T0         = p_spm.tref;    % reference time bin, should be same as ref slice.       
        SPM.xBF.UNITS      = p_spm.time;    % OPTIONS: 'scans'|'secs' for onsets 
        SPM.xBF.Volterra   = 1;             % OPTIONS: 1|2 = order of convolution; 1 = no Volterra
        % Check for presence of FIR variables
        if isfield(p_spm,'hrfwindow') && isfield(p_spm,'hrfbasis')
            SPM.xBF.length     = p_spm.hrfwindow;            % Length of HRF window
            SPM.xBF.order      = p_spm.hrfbasis;             % Order of Basis Set
        end
        if p_spm.design_only == 0
            % Specify data if not design_only
            scans = {};
            scans_files = {};
            for run = p_spm.include_run
                % Check if user wants to use smoothed or unsmoothed files
                % If unsmoothed then check if the fmri token is set to smooth
                % (catched possibility that user never smoothed to begin with)
                % is fmri token is smooth, then take the next step. Otherwise
                % assume user never smoothed and proceed with the fmri token.
                if p_spm.use_smooth == 1
                    epifile = sprintf('%s_r%02d.nii',p_spm.fmri,run);
                elseif p_spm.fmri(1) == 's'
                    epifile = sprintf('%s_r%02d.nii',p_spm.fmri(2:end),run);
                else
                    epifile = sprintf('%s_r%02d.nii',p_spm.fmri,run);
                end
                % Check that files aren't gzipped and stealth unzip if they are.
                epifilepath = fullfile(p_spm.datadir,epifile);
                if exist([epifilepath,'.gz'],'file') == 2    
                    gunzip([epifilepath,'.gz'])
                end    
                scans{end+1} = spm_select('expand', epifilepath);
                scans_files{end+1} =  epifilepath;
            end        
            % Set to char array (SPM expects char array not cell array)
            SPM.xY.P = char(scans);  
            % Add the fullpath to files without selected volumes. This is a 
            % internal hack so that we can pass this to gunzip later on 
            % in case the files are zipped. In this unlikely even that SPM
            % decides to use this Pfiles fieldname, then all hell will break 
            % loose.
            SPM.xY.Pfiles = char(scans_files);
        end
        % Setup Session Structure Array for events, blocks and regressors.
        SPM.Sess.C.name = {};
        SPM.Sess.C.C    = [];
        SPM.Sess.U      = [];       
        for mfield = {'events','blocks','regressors'};
            if ~isempty(p_spm.(mfield{1}))
                for onsfields = p_spm.(mfield{1})
                    if strcmp(mfield{1},'blocks') && p_spm.block_conv == 0
                        spm12w_logger('msg', ['[DEBUG] block_conv ', ...
                            'parameter = 0. Blocks will not be convolved ', ...
                            'with HRF'],'level',p_spm.loglevel)
                        % Generate custom block onset and add tp C
                        block = zeros(sum(p_spm.nvols),1);
                        % Switch from zero indexing to 1 indexing to adjust
                        % for fact that onsets start at zero.
                        block_ons = p_spm.X_onsets.(onsfields{1}).ons + 1;
                        block_dur = p_spm.X_onsets.(onsfields{1}).dur;
                        for i = 1:length(block_ons)
                            block(block_ons(i):block_ons(i)+block_dur(i),1) = 1;
                        end
                        % Fix any underscores as they cause print issues later
                        SPM.Sess.C.name{end+1} = strrep(onsfields{1},'_','-');
                        SPM.Sess.C.C(:,end+1)  = block;                  
                    elseif strcmp(mfield{1},'regressors')
                        % Fix any underscores as they cause print issues later
                        SPM.Sess.C.name{end+1} = strrep(onsfields{1},'_','-');
                        SPM.Sess.C.C(:,end+1)  = p_spm.X_onsets.(onsfields{1}).ons;                        
                    else %events and regular blocks
                        % Fix any underscores as they cause print issues later
                        SPM.Sess.U(end+1).name = {strrep(onsfields{1},'_','-')};  
                        % Set orthogonalization (spm12 now has option in spm struct)
                        SPM.Sess.U(end).orth = p_spm.orth; %0 = no orth.
                        % Get onsets, durations from X_onsets structure
                        SPM.Sess.U(end).ons  = p_spm.X_onsets.(onsfields{1}).ons;       
                        SPM.Sess.U(end).dur  = p_spm.X_onsets.(onsfields{1}).dur;
                        % Set parametrics
                        % NB. P(ii).h can be zero but is usually 1 when using spm
                        % gui. Therefore I'm setting to 1. 
                        if isfield(p_spm.X_onsets.(onsfields{1}),'P')
                            for ii = 1:length(p_spm.X_onsets.(onsfields{1}).P)
                                SPM.Sess.U(end).P(ii).name = p_spm.X_onsets.(onsfields{1}).P(ii).name;
                                SPM.Sess.U(end).P(ii).P    = p_spm.X_onsets.(onsfields{1}).P(ii).P;
                                SPM.Sess.U(end).P(ii).h    = 1;        % order of polynomial expansion
                                SPM.Sess.U(end).P(ii).i    = [1,ii+1]; % sub-indices of U(i).u for plotting (must start at 2)
                            end
                        else
                            SPM.Sess.U(end).P.name = 'none';
                        end                  
                    end
                end
            end
        end
        for mfield = {'outliers','trends','move'}; %depreciated constants
            if isfield(p_spm,['X_',mfield{1}])
                Csize = size(p_spm.(['X_',mfield{1}]),2);
                Cname = {['r-',mfield{1}]};
                SPM.Sess.C.name(end+1:end+Csize) = repmat(Cname,1,Csize);
                SPM.Sess.C.C(:,end+1:end+Csize)  = p_spm.(['X_',mfield{1}]);
            end
        end     
        
    case 'rfx'
        % Setup SPM Structure for 2nd level RFX 
        spm12w_logger('msg',sprintf(['[DEBUG] Setting up SPM structure for ',...
              'type:%s'], args.type),'level',p_spm.loglevel)
        % Unlike with 1st level designs, we're going to only setup part of 
        % SPM structure ourselves and instead allow spm_run_factorial_design.m 
        % to generate the SPM.mat file for 2nd level designs.
        switch p_spm.rfx_type 
            case 'one-sample'
                job.dir = {p_spm.rfxcondir}; % needs to be dir in cell
                job.des.t1.scans = p_spm.rfxconfiles;
                job.cov = struct('c',{},'cname',{},'iCFI',{},'iCC',{});
                job.multi_cov = struct('files',{},'iCFI',{},'iCC',{});
                job.masking.tm.tm_none = 1;
                job.masking.im = p_spm.rfx_im;
                job.masking.em = {''};
                job.globalc.g_omit = 1;
                job.globalm.gmsca.gmsca_no = 1;
                job.globalm.glonorm = 1;
                [~,dirname] = fileparts(p_spm.rfxcondir); % get the conname for msg.
            case 'two-sample'
                % do nothing for now
            case 'anova1'
                job.dir = {p_spm.rfxdir}; % needs to be dir in cell
                job.des.fblock.fac(1).name = 'subject';
                job.des.fblock.fac(1).dept = 0;
                job.des.fblock.fac(1).variance = 0;
                job.des.fblock.fac(1).gmsca = 0;
                job.des.fblock.fac(1).ancova = 0;
                job.des.fblock.fac(2).name = p_spm.rfx_name;
                job.des.fblock.fac(2).dept = 1;
                job.des.fblock.fac(2).variance = 0;
                job.des.fblock.fac(2).gmsca = 0;
                job.des.fblock.fac(2).ancova = 0;
                job.des.fblock.fsuball.specall.scans = p_spm.rfxconfiles;
                job.des.fblock.fsuball.specall.imatrix = p_spm.rfxfacmat;
                job.des.fblock.maininters{1}.fmain.fnum = 1;
                job.des.fblock.maininters{2}.fmain.fnum = 2;
                job.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
                job.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
                job.masking.tm.tm_none = 1;
                job.masking.im = p_spm.rfx_im;
                job.masking.em = {''};
                job.globalc.g_omit = 1;
                job.globalm.gmsca.gmsca_no = 1;
                job.globalm.glonorm = 1;
                [~,dirname] = fileparts(p_spm.rfxdir); % get the oneway name for msg.                                     
            otherwise
                spm12w_logger('msg',sprintf(['[EXCEPTION] Unknown rfx ' ... 
                'analysis type: %s'], p_spm.rfx_type),'level',p_spm.loglevel)
                error('Unknown rfx analysis type: %s', p_spm.rfx_type)        
        end
        % Fill in the spm structure 
        % Note the the spm12/config dir must be in path, this is set in
        % spm12w_glm_rfx.m         
        spm12w_logger('msg',sprintf(['Generating 2nd level rfx model (%s) ',...
                      'on contrast: %s'], p_spm.rfx_type, dirname),...
                      'level',p_spm.loglevel)
        spmloc = spm_run_factorial_design(job);
        % load the generated SPM file
        SPM = load(spmloc.spmmat{1}); % load the generated SPM mat 
                                      % so we can return it
        SPM = SPM.SPM; % remove nested SPM structure  
        
    otherwise
        spm12w_logger('msg',sprintf(['[EXCEPTION] Unknown SPM structure ' ... 
                  'type:%s'], args.type),'level',p_spm.loglevel)
        diary off
        error('Unknown spm structure type: %s', args.type)    
end