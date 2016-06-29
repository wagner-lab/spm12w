function spm12w_prepare(varargin)
% spm12w_prepare(scannerid, sid, rawformat)
%
% Inputs
% ------
% scannerid: Filename of archived raw scanner data to be "prepared" for spm12w.
%            If left unspecified, then spm12w_prepare will open a window
%            allowing the user to select the scannerid to prepare based on the
%            information found in the subid_list.txt or subid_list file 
%            (see spm12w wiki for more information on this file). 
%
% sid:       sid name that the raw scanner data will be renamed to during
%            conversion. This corresponds to the subject identifier to be used
%            for this study (e.g., s01, s02, s03, etc.). 
%
% rawformat: Format of raw data to prepare. Options are 'nifti', 'parrec' 
%            (Philips) or 'dicom' (Siemens). If left unspecified,
%            spm12w_prepare use the format specifed in the third column of the
%            subid_list.txt            
%
% Converts raw philips format par/rec files or dicom data into nifti format
% files. This will also take nifti files exported by the scanner (e.g. Philips)
% and renames them according to spm12w conventions. All renamed data are
% copied to the study's raw directory. 
%
% spm12w_prepare should be run from the study's root directory, where it
% will look for a subid_list.txt file that specifies the correspondance
% between scanner IDs (i.e., the name of the archive of raw data) and subject
% IDs (i.e., the anonymized name of the subject for the current study). In
% addition spm12w_prepare expects the raw data to be stored in
% parrec_scannerID.tar.gz files for par/rec or nifti_scannerID.tar.gz files for
% scanner nifti files or dicom_scannerID.tar.gz files for DICOM files. Zip 
% archives of directories (e.g., parrec_scannerID.gz or dicom_scannerID.zip) 
% are also supported.
% 
% spm12w_prepare workflow:
%   - Unarchive the raw data (par/rec or nifti or dicom)
%   - Convert to nifti (if par/rec or dicom)
%   - Copy to raw directory using the provided sid for the directory name.
%   - Detect anatomy and epi files and rename them from scanner naming
%     conventions to spm12w conventions.
%
% Additional details on par/rec conversion: 
% PAR/REC data from 2011 Phillips was stored data differently than 2012. The 
% order of slices volumes in the header are different. Thus neither dcm2nii nor 
% the parrec2nii command line tool work. R2AGUI appears to do the job. Also note
% that the parrec2nii par/rec converter gets the anat (or epi) transformation 
% wrong so that they are in poor alignment. r2agui appears to work better. 
% -ddw 2012
%
% Examples:
%
% Without arguments, spm12w will prompt the user to select scanner archive
% files for "preparing" and will rename files and convert rawformats based
% on the information in the subid_list.txt file.
%
%   >> spm12w_prepare
% 
% With arguments, spm12w prepare the specified archive files and convert
% according to the specified rawformat. If the subid_list.txt file contains
% a rawformat, the argument provided here will override that format. Filenames
% do not need to include the dicom/parrec/nifti prefix (see example below). 
%
% This case is useful if you have not yet made a subid_list file. 
%
%   >> spm12w_prepare('scannerid',{'01jan16aa',01jan16bb'}, ...
%                     'sid', {'s01','s02'},'rawformat','dicom')
%
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: March, 2013 | Updated: May, 2016
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
arg_defaults = struct('scannerid','', 'sid', '', 'rawformat','');
args = spm12w_args('nargs',0, 'defaults',arg_defaults, 'arguments',varargin);

% Paths
root = pwd;
archdir = fullfile(root, 'arch');
rawdir = fullfile(root, 'raw');

if isempty(args.scannerid) && isempty(args.sid) && isempty(args.rawformat)
    [scannerlist, sids, rawformats, excludeseries] = spm12w_prepare_csvparser('archdir',...
                                                archdir, 'rootdir', root, ...
                                                'scannerid', args.scannerid, ...
                                                'rawformat', args.rawformat);
else
    % Check for cell in case user provided string
    if ~iscell(args.scannerid)
        args.scannerid = cellstr(args.scannerid);
    end
    if ~iscell(args.sid)
        args.sid = cellstr(args.sid);
    end
    if ~iscell(args.rawformat)
        args.rawformat = cellstr(args.rawformat);
    end 
    scannerlist = args.scannerid;
    sids = args.sid;
    rawformats = args.rawformat;
end

% set R2AGUI options regardless of rawformat
options.subaan = 1;
options.usealtfolder = 0;
options.altfolder = '~/';
options.prefix= [];
options.angulation= 1;
options.rescale= 1;
options.usefullprefix= 0;
options.outputformat= 1;
options.dim= 4;
options.dti_revertb0= 0;

% Prepare either parrec or nifti or dicom files
for i = 1:length(scannerlist)
    % Assign directories and make a tempdir for extracting subjects.
    datadir = fullfile(archdir, rawformats{i});
    subjpath = fullfile(tempname); 
    % Find file and announce intentions
    archfile = dir(fullfile(datadir, ...
                   sprintf('%s_%s.*', rawformats{i},scannerlist{i})));
    % If archfile is empty than either sid is missing from arch or there's a
    % typo in the arch filename or the subid_list. Either way give an error.
    if isempty(archfile)
        spm12w_logger('msg',sprintf(['[EXCEPTION] Cannot find file ',...
                      '%s_%s in %s... Are you sure it exists?'],...
                      rawformats{i}, scannerlist{i},datadir))
        error('Cannot find file %s_%s... Are you sure it exists?', ...
              rawformats{i},scannerlist{i})
    else
        spm12w_logger('msg',sprintf('Subject: %s (%s)',scannerlist{i}, sids{i}))
        spm12w_logger('msg',sprintf(['Sid:%s | Extracting %s data from ',...
                      'archive this may take awhile...'], sids{i},rawformats{i}))
    end    
    % Figure out file extension. Because we can encoutner .tar.gz we can't
    % use fileparts and have to search the string. Next, extract to tmp folder
    % defined above. 
    if ~isempty(strfind(archfile(1).name, '.tar.gz'))
        mkdir(subjpath)
        untar(fullfile(datadir, ...
              sprintf('%s_%s.tar.gz',rawformats{i},scannerlist{i})), subjpath);
    elseif ~isempty(strfind(archfile(1).name, '.gz'))
        mkdir(subjpath)
        gunzip(fullfile(datadir, ...
               sprintf('%s_%s.gz',rawformats{i},scannerlist{i})), subjpath);
    elseif ~isempty(strfind(archfile(1).name, '.zip'))
        mkdir(subjpath)
        unzip(fullfile(datadir, ...
              sprintf('%s_%s.zip',rawformats{i},scannerlist{i})),subjpath);
    else
        spm12w_logger('msg',sprintf(['[EXCEPTION] Sid:%s | Unreocgnized ',...
                      'archive extension (not zip, gz or tar.gz)...'],sids{i}))
        error('Unrecognized archive extension. Aborting...');
    end
    if strcmpi(rawformats{i}, 'nifti') 
        % Find the raw NIFTI files for this subject
        files = dir(fullfile(subjpath,scannerlist{i},'*.nii'));
        json=0; % no json metadata with scanner nifti files unless we decide to extract some.
    elseif strcmpi(rawformats{i}, 'parrec')       
        % Find the raw PAR files for this subject
        files = dir(fullfile(subjpath,scannerlist{i},'*.PAR'));
        filelist = {files.name};
        options.pathpar = fullfile(subjpath,scannerlist{i},filesep);
        % Convert from parrec if rawformat is parrec
        spm12w_logger('msg','Converting format: parrec...')
        convert_r2a(filelist,options);
        json = 0; %at the moment no json files, check convert_r2a if we can grab some metadata.
    elseif strcmpi(rawformats{i},'dicom')
        spm12w_logger('msg','Converting format: dicom...')
        % Set preferences for dicm2nii. 
        setpref('dicm2nii_gui_para', 'save_json', true); %save headers as json
        setpref('dicm2nii_gui_para', 'use_parfor', false); %disable dicm2nii parfor
        setpref('dicm2nii_gui_para', 'save_patientName', false); %save id       
        dicm2nii(subjpath, subjpath, 'nii');       
        % Find the converted NIFTI files for this subject
        files = dir(fullfile(subjpath,'*.nii'));
        json = 1; %json files exist.
    else
        spm12w_logger('msg',sprintf(['[EXCEPTION] Sid:%s | Unknown raw ',...
                      'format: %s...'],sids{i}, rawformats{i}));
        error('Unknown rawformat... Aborting...')
    end
    
    % Trim the excluded series, fix the series that we keep and trim 
    % the small stuff.
    spm12w_logger('msg',sprintf(['Sid:%s | Pruning unecessary files',...
                                 ' from conversion'],sids{i}));
    count = 1;
    for ii = 1:length(files)
        % Prune scout and other uncessary scans below 10MB.
        % But be mindful of fieldmaps. 
        % At some point we should harmonize these hacks across
        % dicom/parrec/nifti and the various sites. 
        if isempty(strfind(files(ii).name,sprintf('_s%03d.nii',str2num(excludeseries{i}))))
            if files(ii).bytes > 10000000 || ~isempty(strfind(files(ii).name,'fieldmap'))
                % Check for the series that we keep if there's excludes
                if ~isempty(excludeseries{i}) && ~isempty(regexp(files(ii).name,'_s0\d\d.nii', 'once'))
                    % if true rename file and it json and add to the filelist
                    spm12w_logger('msg',sprintf('[DEBUG] Sid:%s | Adjusting name of unexcluded series: %s',...
                          sids{i},files(ii).name));
                    newname = sprintf('%s.nii',files(ii).name(1:regexp(files(ii).name,'_s0\d\d.nii')-1));
                    spm12w_logger('msg',sprintf('[DEBUG] Sid:%s | Renaming %s to %s...',...
                          sids{i},files(ii).name, newname));
                    movefile(fullfile(subjpath,files(ii).name),fullfile(subjpath,newname))
                    [~,oldnamejson]=fileparts(files(ii).name);
                    oldnamejson = [oldnamejson,'.json'];
                    newnamejson = sprintf('%s.json',files(ii).name(1:regexp(files(ii).name,'_s0\d\d.nii')-1));
                    spm12w_logger('msg',sprintf('[DEBUG] Sid:%s | Renaming %s to %s...',...
                          sids{i},oldnamejson, newnamejson));
                    movefile(fullfile(subjpath,oldnamejson),fullfile(subjpath,newnamejson))
                    filelist{count} = newname;
                    count = count + 1;
                else
                    filelist{count} = files(ii).name;
                    count = count + 1;
                end
            else
                spm12w_logger('msg',sprintf('[DEBUG] Sid:%s | Skipping file: %s',...
                              sids{i},files(ii).name));
            end
        else
            spm12w_logger('msg',sprintf('[DEBUG] Sid:%s | Skipping file: %s',...
                          sids{i},files(ii).name));
        end
    end   
       
    % Make various subject dirs
    rawsid = fullfile(rawdir,sids{i});
    if exist(rawsid,'dir')
        spm12w_logger('msg',sprintf(['[DEBUG] Sid:%s | Previous Subject ',...
                      'directory exists...'], sids{i}))
        spm12w_logger('msg',sprintf(['[DEBUG] Sid:%s | Previous directory ',...
                      'exists, deleting and preparing...'], sids{i}))        
        % Delete prior subject directory
        rmdir(rawsid,'s')
    end
    mkdir(rawsid);
    
    % Move files
    for i_flist = 1:size(filelist,2)
        [~,file_noext] = fileparts(filelist{i_flist});
        if strcmpi(rawformats{i}, 'parrec')     
            moveme = fullfile(subjpath,scannerlist{i},file_noext);
            movefile([moveme,filesep,'*.nii'],rawsid);       
        elseif strcmpi(rawformats{i}, 'nifti')
            moveme = fullfile(subjpath,scannerlist{i},[file_noext,'.nii']);
            movefile(moveme,rawsid);   
        elseif strcmpi(rawformats{i}, 'dicom')
            moveme = fullfile(subjpath,[file_noext,'.nii']);
            movefile(moveme,rawsid);  
            moveme = fullfile(subjpath,[file_noext,'.json']);
            movefile(moveme,rawsid);              
        end
        % Check for bval and vec file (OSU data)
        if exist(fullfile(subjpath,'DTI_64.bval'),'file')
            movefile(fullfile(subjpath,'DTI_64.bval'),...
                     fullfile(rawsid,'dti.bval'));      
        end
        if exist(fullfile(subjpath,'DTI_64.bvec'),'file')
            movefile(fullfile(subjpath,'DTI_64.bvec'),...
                     fullfile(rawsid,'dti.bvec'));      
        end
    end   
    
    % Rename files
    if strcmpi(rawformats{i}, 'parrec') 
        % Anatomical renaming
        mrifile = dir(fullfile(rawsid,'*T1TFE*.nii'));
        if ~isempty(mrifile)
            spm12w_logger('msg',sprintf(['[DEBUG] Sid:%s | Renaming %s to ',...
                          'anat.nii and compressing...'], sids{i}, mrifile(1).name))
            movefile(fullfile(rawsid,mrifile(1).name),fullfile(rawsid,'anat.nii'));
            gzip(fullfile(rawsid,'anat.nii'));
            delete(fullfile(rawsid,'anat.nii'));
            if size(mrifile,1) == 2
                movefile(fullfile(rawsid,mrifile(2).name),fullfile(rawsid,'anat_HalfAndHalf.nii'));
                gzip(fullfile(rawsid,'anat_HalfAndHalf.nii'));
                delete(fullfile(rawsid,'anat_HalfAndHalf.nii'));
            end
        else
            if exist(fullfile(datadir,[scannerlist{i},'_anat'],'anat.nii.gz'),'file')
                spm12w_logger('msg',sprintf(['[DEBUG] Sid:%s | Copying ',...
                              'anonymized anatomical...'], sids{i}))
                copyfile(fullfile(datadir,[scannerlist{i},'_anat'],'anat.nii.gz'), ...
                    rawsid)
            end
        end
        % DTI moving and renaming
        dtifiles = dir(fullfile(rawsid,'*DwiSE*.nii'));
        if ~isempty(dtifiles)
            spm12w_logger('msg',sprintf(['[DEBUG] Sid:%s | Renaming and ',...
                          'compressing dti files...'], sids{i}));
            dtirun = 1;
            lastgradnum = 0;
            for x = 1:length(dtifiles)
                %-Figure out the appropriate DTI name (2 runs plus the g
                %-numbers)
                gradnum = dtifiles(x).name;
                gradnum = gradnum(strfind(gradnum,'-g')+2:end-4);
                gradnum = str2num(gradnum);
                if gradnum < lastgradnum
                    dtirun = dtirun + 1;
                end
                lastgradnum = gradnum;
                dtiname = sprintf('dti%d_g%02d.nii',dtirun,gradnum);
                movefile(fullfile(rawsid,dtifiles(x).name),fullfile(rawsid,dtiname));
                gzip(fullfile(rawsid,dtiname));
                delete(fullfile(rawsid,dtiname));
            end
        end
        % epi file renaming
        boldfiles = dir(fullfile(rawsid,'*FEEPI*.nii'));
        epi_count = 1;
        for x = 1:length(boldfiles)
            boldname = sprintf('epi_r%.2d.nii',epi_count);
            epi_count = epi_count + 1;
            spm12w_logger('msg',sprintf(['[DEBUG] Sid:%s | Renaming %s to ',...
                          '%s and compressing...'],sids{i},...
                          boldfiles(x).name,boldname));
            movefile(fullfile(rawsid,boldfiles(x).name),fullfile(rawsid,boldname));
            gzip(fullfile(rawsid,boldname));
            delete(fullfile(rawsid,boldname));
        end
    elseif strcmpi(rawformats{i}, 'nifti') 
        niifiles = dir(fullfile(rawsid,'*.nii'));
        epi_count = 1;
        anat_count = 1;
        for x = 1:length(niifiles)
            if niifiles(x).bytes > 30000000 %anatomicals are under 30 megs
                boldname = sprintf('epi_r%.2d.nii', epi_count);
                epi_count = epi_count + 1;
                spm12w_logger('msg',sprintf(['[DEBUG] Sid:%s | Renaming %s ',...
                              'to %s and compressing...'],sids{i}, ...
                              niifiles(x).name,boldname));
                movefile(fullfile(rawsid,niifiles(x).name),fullfile(rawsid,boldname));
                gzip(fullfile(rawsid,boldname));
                delete(fullfile(rawsid,boldname));
            end    
            % Anatomical renaming
            if niifiles(x).bytes > 28800000 && niifiles(x).bytes < 29000000
                if anat_count == 1
                    spm12w_logger('msg',sprintf(['[DEBUG] Sid:%s | Renaming %s to anat.nii ',...
                                  'and compressing'],sids{i},niifiles(x).name));
                    movefile(fullfile(rawsid,niifiles(x).name),fullfile(rawsid,'anat.nii'));
                    gzip(fullfile(rawsid,'anat.nii'));
                    delete(fullfile(rawsid,'anat.nii'));
                    anat_count = anat_count + 1;
                elseif anat_count == 2
                    spm12w_logger('msg',sprintf(['[DEBUG] Sid:%s | Renaming %s to anat_2.nii ',...
                                  'and compressing'], sids{i},niifiles(x).name));
                    movefile(fullfile(rawsid,niifiles(x).name),fullfile(rawsid,'anat_2.nii'));
                    gzip(fullfile(rawsid,'anat_2.nii'));
                    delete(fullfile(rawsid,'anat_2.nii'));
                    anat_count = anat_count + 1;
                end   
            end
        end
    elseif strcmpi(rawformats{1}, 'dicom')
        niifiles = dir(fullfile(rawsid,'*.nii'));
        epi_count = 1;
        rest_count = 1;
        % Naming conventions are closer to spm12w on OSU dicoms making this easy
        for niifile = {niifiles.name}
            switch logical(true)          
                case strcmpi(niifile{1},sprintf('run%d.nii',epi_count))
                    newname = sprintf('epi_r%.2d.nii', epi_count);
                    epi_count = epi_count + 1;
                case strcmpi(niifile{1},sprintf('rest%d.nii',rest_count))
                    newname = sprintf('rest_r%.2d.nii', rest_count);
                    rest_count = rest_count + 1;
                case strcmpi(niifile{1},'MPRAGE.nii')
                    newname = 'anat.nii';
                case strcmpi(niifile{1},'DTI_64.nii')
                    newname = 'dti.nii';                 
                case strcmpi(niifile{1},'fieldmap.nii')
                    newname = 'epi_fieldmap.nii';
                case strcmpi(niifile{1},'fieldmap_phase.nii')
                    newname = 'epi_fieldmap_phase.nii';        
                otherwise
                    % cleanup tmp dir
                    rmdir(subjpath,'s');
                    spm12w_logger('msg',sprintf(['[EXCEPTION] Sid:%s | ',...
                                   'Unknown filetype: %s...'],sids{i},niifile{1}));
                    error('unknown filetype. Aborting...')
            end
            spm12w_logger('msg', sprintf(['[DEBUG] Sid:%s | Renaming %s to ',...
                          '%s and compressing...'],sids{i},niifile{1},newname));
            movefile(fullfile(rawsid,niifile{1}),fullfile(rawsid,newname));
            gzip(fullfile(rawsid,newname));
            delete(fullfile(rawsid,newname));     
            % now for json
            jsonfile = strrep(niifile{1},'.nii','.json');
            newname = strrep(newname,'.nii','.json');
            spm12w_logger('msg', sprintf(['[DEBUG] Sid:%s | Renaming %s to ',...
                          '%s...'],sids{i},jsonfile,newname));
            movefile(fullfile(rawsid,jsonfile),fullfile(rawsid,newname));           
        end
    end
    % Delete the extracted data
    rmdir(subjpath,'s');
end

% Print completed
spm12w_logger('msg', sprintf(['Conversion to nifti on %d subjects ',...
              'completed...'], length(scannerlist)));
cd(root)