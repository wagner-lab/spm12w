function spm12w_dirsetup(varargin)
% spm12w_dirsetup(dirtype, params)
%
% Input
% -----
% dirtype : Type of directory to setup (e.g., 'prep_clean', 'prep_setup', 
%           'glm_clean', 'glm_setup'). (Default='prep_clean').
%
% params  : A structure of parameters (e.g., p, glm, etc.). 
%
% Takes a parameters structure and a dirtype and prepares the directories
% for the desired stage of analysis or preprocessing. Dirtype can be:
%
%   'prep_clean' : spm12w_dirsetup will check for a pre-existing directory and 
%                  will preseve log files by moving them to an archive directory
%                  with a date and time stamp. All other files and directories
%                  within the preprocessing data directory will be deleted.
%                  Any prior preprocessing directories will be deleted if they
%                  exist. If volume outliers were created, these will also be
%                  deleted as they should be regenerated whenever the 
%                  preprocessing parameters change. NB: if people request it, 
%                  this could be changed in future versions, in which case
%                  outlier files will NOT be deleted).
%                  If a prior preprocessing directory exists, its preprocess 
%                  log, mat and pdf files will be saved in an archive directory
%                  in the preprocesisng dir with the date and time prepended. 
%                  Finally, this function will check for standard spm12w 
%                  directories and create them if they do not exist (i.e., 
%                  auxil, notes, qa).
%
%   'prep_setup' : spm12w_dirsetup copy raw nifti data to the p.datadir for 
%                  preprocessing, renaming them appropriately.
%
%   'glm_clean'  : spm12w_dirsetup will check for a pre-existing directory and 
%                  will preseve log files by moving them to an archive directory
%                  with a date and time stamp. All other files and directories
%                  within the glm analysis directory will be deleted.
%
%   'con_clean'  : spm12w_dirsetup will check the pre-existing glm directory and 
%                  remove all prior con and spmT files (to prevent outdated
%                  con files from lingering around). 
%
%   'rfx_clean'  : spm12w_dirsetup will check the pre-existing rfx directory and 
%                  remove all prior files associated with the rfx analysis.
%                  The prior spmT and SPM.mat files will be preserved by 
%                  moving them to an archive directory with a date and time
%                  stamp.
%
% Examples:
%
%   >> spm12w_dirsetup('dirtype', 'prep_clean', 'params', p)
%   >> spm12w_dirsetup('dirtype', 'glm_clean', 'params', glm)
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: November, 2014 | Updated: April, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('dirtype','prep', 'params','');
args = spm12w_args('nargs',4, 'defaults', args_defaults, 'arguments', varargin);

% Assign params to p for cleaner code. 
dirp = args.params;

% Setup directories
spm12w_logger('msg',dirp.niceline, 'level',dirp.loglevel)

switch args.dirtype
    case 'prep_clean'
        spm12w_logger('msg',sprintf('Checking prep directory for subject: %s', ...
                      dirp.sid),'level',dirp.loglevel)
        % Check that studyroot contains spm12w directories (i.e., arch, raw and
        % scripts should already exist or else it's not possible to run
        % preprocessing). We don't need to check for the datadir since
        % we will create it.
        spm12w_logger('msg','Checking study root paths...','level',dirp.loglevel)
        dirlist = {dirp.anadir, dirp.auxdir, dirp.notesdir, dirp.qadir};
        for i = 1:length(dirlist)
           if ~isdir(dirlist{i})
                spm12w_logger('msg',sprintf('[DEBUG] %s not found. Creating...',...
                              dirlist{i}),'level',dirp.loglevel)
                mkdir(dirlist{i})
           end 
        end
        % Check for prior preprocessing and make if missing.
        if ~exist(dirp.datadir,'dir')
            spm12w_logger('msg',sprintf('[DEBUG] Creating directory: %s ',...
              dirp.datadir),'level',dirp.loglevel)
            mkdir(dirp.datadir)
        else
            spm12w_logger('msg',['[WARNING] Prior prep directory found. ', ...
                          'Directory will be deleted...'],'level',dirp.loglevel)
            spm12w_logger('msg',sprintf(['[DEBUG] Deleting prior ',...
                          'prep directory at: %s'], dirp.datadir), ...
                          'level',dirp.loglevel) 
            flist = dir(dirp.datadir);
            flist = flist(3:end); % remove . and ..
            timestamp = datestr(now, 'dd-mm-yyyy_HH-MM'); %timestamp for renaming
            for i = 1:length(flist) 
                diary off % in case a logger is still open (diary stops rmdir)
                if flist(i).isdir
                    if strcmp(flist(i).name,'archive')
                        spm12w_logger('msg','[DEBUG] Preserving archive directory', ...
                                      'level', dirp.loglevel) % Placeholder
                    else
                        spm12w_logger('msg',sprintf('[DEBUG] Removing directory: %s',...
                                      flist(i).name),'level',dirp.loglevel)
                        rmdir(fullfile(dirp.datadir,flist(i).name))
                    end
                else
                    % Keep certain files
                    [~,fname,~] = fileparts(dirp.preplog);
                    keepfiles = {[fname,'.log'],[fname,'.pdf'],[fname,'.mat']};               
                    if ismember(flist(i).name,keepfiles)
                        if ~exist(dirp.preparch,'dir')
                            mkdir(dirp.preparch)
                        end
                        infile = fullfile(dirp.datadir,flist(i).name);
                        outfile = fullfile(dirp.preparch, ... 
                                  sprintf('%s_%s',timestamp,flist(i).name));
                        spm12w_logger('msg',sprintf(['[DEBUG] Prior logfile ', ...
                                     'found. Moving %s to %s'],infile, outfile), ...
                                     'level',dirp.loglevel)
                        movefile(infile,outfile);
                    else
                        spm12w_logger('msg',sprintf('[DEBUG] Removing file: %s', ...
                                      flist(i).name),'level',dirp.loglevel)
                        delete(fullfile(dirp.datadir,flist(i).name))
                    end            
                end
            end
        end
        
    case 'prep_setup'
        spm12w_logger('msg',sprintf('Setting up prep directory for subject: %s', ...
              dirp.sid),'level',dirp.loglevel)
        % Copy raw files over to prep dir
        mrifile = sprintf('%s.nii.gz', dirp.mri);
        if ~exist(fullfile(dirp.rawdir,mrifile),'file')
            msg = sprintf(['[WARNING] Missing %s file for subject: %s... ' ...
                          'proceeding...'],mrifile,dirp.sid);
            spm12w_logger('msg',msg,'level',dirp.loglevel)
        else
            spm12w_logger('msg',sprintf('Copying and extracting file: %s...',mrifile),...
                          'level',dirp.loglevel)
            gunzip(fullfile(dirp.rawdir,mrifile),dirp.datadir)
        end

        % Copy epi files over to prep dir
        for i = 1:dirp.nses
            niifile = sprintf('%s_r%02d.nii.gz', dirp.fmri, i);
            if ~exist(fullfile(dirp.rawdir,niifile),'file')
                spm12w_logger('msg',sprintf(['[EXCEPTION] Missing %s file for ' ... 
                          'subject: %s'], niifile, dirp.sid),'level',dirp.loglevel)
                diary off
                error('Missing %s file for subject: %s', niifile, dirp.sid)
            else
                spm12w_logger('msg',sprintf('Copying and extracting file: %s...',...
                              niifile),'level',dirp.loglevel)
                gunzip(fullfile(dirp.rawdir,niifile),dirp.datadir)   
            end
        end

    case 'glm_clean'
        spm12w_logger('msg',sprintf('Checking glm directory for subject: %s', ...
                      dirp.sid),'level',dirp.loglevel)
        if ~exist(dirp.glmdir,'dir')
            spm12w_logger('msg',sprintf('[DEBUG] Creating directory: %s ',...
              dirp.glmdir),'level',dirp.loglevel)
            mkdir(dirp.glmdir)
        else
            spm12w_logger('msg',['[WARNING] Prior glm directory found. ', ...
                          'Directory will be cleaned...'],'level',dirp.loglevel)
            spm12w_logger('msg',sprintf(['[DEBUG] Cleaning prior ',...
                          'glm directory at: %s'], dirp.glmdir), ...
                          'level',dirp.loglevel)  
            flist = dir(dirp.glmdir);
            flist = flist(3:end); % remove . and ..
            timestamp = datestr(now, 'dd-mm-yyyy_HH-MM'); %timestamp for renaming
            for i = 1:length(flist) 
                diary off % in case a logger is still open (diary stops rmdir)
                if flist(i).isdir
                    if strcmp(flist(i).name,'archive')
                        spm12w_logger('msg','[DEBUG] Preserving archive directory', ...
                                      'level', dirp.loglevel) % Placeholder
                    else
                        spm12w_logger('msg',sprintf('[DEBUG] Removing directory: %s',...
                                      flist(i).name),'level',dirp.loglevel)
                        rmdir(fullfile(dirp.glmdir,flist(i).name))
                    end
                else
                    % Keep log, pdf and mat files
                    [~,fname,~] = fileparts(dirp.glmlog);
                    keepfiles = {[fname,'.log'],[fname,'.pdf'],[fname,'.mat']};               
                    if ismember(flist(i).name,keepfiles)
                        if ~exist(dirp.glmarch,'dir')
                            mkdir(dirp.glmarch)
                        end
                        infile = fullfile(dirp.glmdir,flist(i).name);
                        outfile = fullfile(dirp.glmarch, ... 
                                  sprintf('%s_%s',timestamp,flist(i).name));
                        spm12w_logger('msg',['[DEBUG] Prior logfile ', ...
                                     'found.'], 'level',dirp.loglevel)
                        spm12w_logger('msg',sprintf('[DEBUG] Moving %s',...
                                     infile),'level',dirp.loglevel)
                        spm12w_logger('msg',sprintf('[DEBUG] To %s',...
                                     outfile),'level',dirp.loglevel)
                        movefile(infile,outfile);
                    else
                        spm12w_logger('msg',sprintf('[DEBUG] Removing file: %s', ...
                                      flist(i).name),'level',dirp.loglevel)
                        delete(fullfile(dirp.glmdir,flist(i).name))
                    end            
                end
            end
        end
       
    case 'con_clean'
        spm12w_logger('msg',sprintf('Checking glm directory for subject: %s', ...
                      dirp.sid),'level',dirp.loglevel)
        if ~exist(dirp.glmdir,'dir')
            spm12w_logger('msg',sprintf('[WARNING] Cannot find prior glm directory: %s ',...
              dirp.glmdir),'level',dirp.loglevel)
            error('Cannot find prior glm directory at: %s', dirp.glmdir) 
        else
            for fprefix = {'con_0*.nii','spmT_0*.nii','spmF_*.nii','ess_0*.nii'}; 
                flist = ls(fullfile(dirp.glmdir,fprefix{1}));
                if ~isempty(flist)
                    spm12w_logger('msg',sprintf(['[WARNING] Prior contrast files (%s) found. ', ...
                          'Previous contrast files will be deleted...'],fprefix{1}),'level',dirp.loglevel) 
                    for f_i = 1:size(flist,1)
                        spm12w_logger('msg',sprintf('[DEBUG] Removing file: %s', ...
                                      deblank(flist(f_i,:))),'level',dirp.loglevel)
                        delete(fullfile(dirp.glmdir,deblank(flist(f_i,:))))
                    end
                end     
            end
        end
   case 'rfx_clean'
    % For rfx dirs, we need to make them on the fly based on dirp.rfx_conds var
    for rfxcondir = dirp.rfx_conds
        rfxdir = fullfile(dirp.rfxdir,rfxcondir{1});
        spm12w_logger('msg',sprintf('Checking rfx directory: %s', ...
                      rfxdir),'level',dirp.loglevel)
        if ~exist(rfxdir,'dir')
            spm12w_logger('msg',sprintf('[DEBUG] Creating directory: %s ',...
              rfxdir),'level',dirp.loglevel)
            mkdir(rfxdir)
        else
            spm12w_logger('msg',['[WARNING] Prior rfx directory found. ', ...
                          'Directory will be cleaned...'],'level',dirp.loglevel)
            spm12w_logger('msg',sprintf(['[DEBUG] Cleaning prior ',...
                          'rfx directory at: %s'], rfxdir), ...
                          'level',dirp.loglevel)  
            rfxarch = fullfile(rfxdir,dirp.rfxarchtok); % set the arch dir name
            timestamp = datestr(now, 'dd-mm-yyyy_HH-MM'); %timestamp for renaming
            % Keep spmT files
            %(todo: keep spm_fs for anova once we add
            % anova to supported rfx models)
            spmtlist = ls(fullfile(rfxdir,'spmT_*.nii'));   
            for f_i = 1:size(spmtlist,1)
                % Create the archive directory if it doesn't exist
                if ~exist(rfxarch,'dir')
                    mkdir(rfxarch)
                end
                infile = fullfile(rfxdir,deblank(spmtlist(f_i,:)));
                outfile = fullfile(rfxarch, ... 
                          sprintf('%s_%s',timestamp,deblank(spmtlist(f_i,:))));
                spm12w_logger('msg',['[DEBUG] Prior analysis ', ...
                             'spmT file found.'], 'level',dirp.loglevel)
                spm12w_logger('msg',sprintf('[DEBUG] Moving %s',...
                             infile),'level',dirp.loglevel)
                spm12w_logger('msg',sprintf('[DEBUG] To %s',...
                             outfile),'level',dirp.loglevel)
                movefile(infile,outfile)
            end
            % Check all other files (keep .mat files, delete others)
            flist = dir(rfxdir);
            flist = flist(3:end); % remove . and ..
            for i = 1:length(flist) 
                if flist(i).isdir
                    if strcmp(flist(i).name,'archive')
                        spm12w_logger('msg','[DEBUG] Preserving archive directory', ...
                                      'level', dirp.loglevel) % Placeholder
                    else
                        spm12w_logger('msg',sprintf('[DEBUG] Removing directory: %s',...
                                      flist(i).name),'level',dirp.loglevel)
                        rmdir(fullfile(rfxdir,flist(i).name))
                    end
                else    
                    % Keep mat files
                    keepfiles = {'SPM.mat',sprintf('%s.mat',dirp.rfx_name)};               
                    if ismember(flist(i).name,keepfiles)
                        % Create the archive directory if it doesn't exist
                        if ~exist(rfxarch,'dir')
                            mkdir(rfxarch)
                        end
                        infile = fullfile(rfxdir,flist(i).name);
                        outfile = fullfile(rfxarch, ... 
                                  sprintf('%s_%s',timestamp,flist(i).name));
                        spm12w_logger('msg',['[DEBUG] Prior analysis ', ...
                                     'mat file found.'], 'level',dirp.loglevel)
                        spm12w_logger('msg',sprintf('[DEBUG] Moving %s',...
                                     infile),'level',dirp.loglevel)
                        spm12w_logger('msg',sprintf('[DEBUG] To %s',...
                                     outfile),'level',dirp.loglevel)
                        movefile(infile,outfile);
                    else
                        spm12w_logger('msg',sprintf('[DEBUG] Removing file: %s', ...
                                      flist(i).name),'level',dirp.loglevel)
                        delete(fullfile(rfxdir,flist(i).name))
                    end            
                end
            end
        end          
    end
        
    otherwise
        spm12w_logger('msg',sprintf('[EXCEPTION] Unknown dirtype: %s', ... 
                      args.dirtype),'level',dirp.loglevel)
        diary off
        error('Unknown dirtype: %s', args.dirtype)    
end

% Print completed message (tailor msg if rfx clean).
if strcmp(args.dirtype,'rfx_clean')
    dirmsg = sprintf(['[DEBUG] Finished setting up directories ',...
                 'for dirtype: %s'], args.dirtype);
else
    dirmsg = sprintf(['[DEBUG] Finished setting up directories ',...
                 'for dirtype: %s and subject: %s'], args.dirtype, dirp.sid);
end
     spm12w_logger('msg',dirmsg,'level',dirp.loglevel)  