function spm12w_prepare(varargin)
% spm12w_prepare(scannerid, rawformat)
%
% Inputs
% ------
% scannerid: Filename of archived raw scanner data to be "prepared" for spm12w.
%            If left unspecified, then spm12w_prepare will open a window
%            allowing the user to select the scannerid to prepare based on the
%            information found in the subid_list.txt file (see spm12w wiki
%            for more information on this file). 
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
%   >> spm12w_prepare('scannerid',{'01jan16aa',01jan16bb'}, ...
%                     'sid', 's01','s02'},'rawformat','dicom')
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
sidfile = fullfile(archdir, 'subid_list.txt');
rawdir = fullfile(root, 'raw');

if isempty(args.scannerid) && isempty(args.sid) && isempty(args.rawformat)
    % Load the subid_list.txt file

    % Check that we are in the study's root directory and that files exists
    if ~exist(sidfile, 'file') 
        error(['I can''t seem to find the subid_list.txt file in the arch '...
            'directory. Are you sure it exists and that %s is the study ' ...
            'root directory?'], root)
    end

    % Load subjects and sids from subid_list.txt file
    filename = sidfile;
    formatSpec = '%s%s%[^\n\r]';
    fid = fopen(filename, 'r');
    subidArray = textscan(fid, formatSpec, 'Delimiter', ',');
    fclose(fid);

    % Set subjects, sids and rawformat to approriate entries in subidArray.
    % transpose to get cell array size to match args.scannerid & args.rawformat.
    scannerlist = deblank(subidArray{1})';
    sids = deblank(subidArray{2})';
    rawformats = deblank(subidArray{3})';

    % Allow user to select the sids they want
    if isempty(args.scannerid)
        scanneridx = listdlg('PromptString', 'Select archive(s):', ...
                              'SelectionMode','multiple','ListString',scannerlist);   
    else
        scanneridx = find(ismember(scannerlist, args.scannerid));    
    end

    % reassign cell arrays according to the idx found above.
    scannerlist = scannerlist(scanneridx);   
    sids = sids(scanneridx);     
    rawformats = rawformats(scanneridx); 
    % If user supplied rawformat arg, then replace rawformat harvested from 
    % subid_list.txt with user supplied rawformat
    if ~isempty(args.rawformat)
        rawformats(:) = {args.rawformat};    
    end
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
    fprintf(['Subject: %s (%s) | Extracting %s data from archive, ', ...
             'this may take awhile...\n'], scannerlist{i}, sids{i}, ...
             rawformats{i});
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
        error(['Unrecognized archive extension (not zip, gz or tar.gz). ', ...
               'Aborting...']);
    end
    if strcmpi(rawformats{i}, 'nifti') 
        % Find the raw NIFTI files for this subject
        files = dir(fullfile(subjpath,scannerlist{i},'*.nii'));
    elseif strcmpi(rawformats{i}, 'parrec')       
        % Find the raw PAR files for this subject
        files = dir(fullfile(subjpath,scannerlist{i},'*.PAR'));
        filelist = {files.name};
        options.pathpar = fullfile(subjpath,scannerlist{i},filesep);
        % Convert from parrec if rawformat is parrec
        fprintf('Converting format: parrec...\n');
        convert_r2a(filelist,options);
        fprintf('Done.\n');
    elseif strcmpi(rawformats{i},'dicom')
        fprintf('Converting format: dicom...\n');
        % Disable parallel looping in xiangrui's script?
        % ps = parallel.Settings;
        % ps.Pool.AutoCreate = false;
        dicm2nii(subjpath, subjpath, 'nii');
        fprintf('Done.\n');
        % Find the converted NIFTI files for this subject
        files = dir(fullfile(subjpath,'*.nii'));
    else
        error('Unknown rawformat: %s. Aborting...', rawformats{i})
    end
    
    % Trim the small stuff
    count = 1;
    for ii = 1:length(files)
        % Prune scout and other uncessary scans below 10MB.
        % But be mindful of fieldmaps. 
        % At some point we should harmonize these hacks across
        % dicom/parrec/nifti and the various sites. 
        if files(ii).bytes > 10000000 && isempty(strfind(files(ii).name,'fieldmap'))
            filelist{count} = files(ii).name;
            count = count + 1;
        elseif ~isempty(strfind(files(ii).name,'fieldmap'))
            filelist{count} = files(ii).name;
            count = count + 1;
        end
    end   
       
    % Make various subject dirs
    rawsid = fullfile(rawdir,sids{i});
    if exist(rawsid,'dir')
        fprintf(['Previous Subject directory exists... deleting the ', ...
                  'directory and preparing from raw.\n'])
        % Delete prior subject directory
        fprintf('Deleting prior subject directory...\n')
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
            fprintf('Renaming %s to %s\n',mrifile(1).name,'anat.nii and compressing...');
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
                fprintf('Copying anonymized anatomical...\n')
                copyfile(fullfile(datadir,[scannerlist{i},'_anat'],'anat.nii.gz'), ...
                    rawsid)
            end
        end
        % DTI moving and renaming
        dtifiles = dir(fullfile(rawsid,'*DwiSE*.nii'));
        if ~isempty(dtifiles)
            fprintf('Renaming and compressing dti files...\n');
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
            fprintf('Renaming %s to %s and compressing...\n',boldfiles(x).name,boldname);
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
                fprintf('Renaming %s to %s and compressing...\n',niifiles(x).name,boldname);
                movefile(fullfile(rawsid,niifiles(x).name),fullfile(rawsid,boldname));
                gzip(fullfile(rawsid,boldname));
                delete(fullfile(rawsid,boldname));
            end    
            % Anatomical renaming
            if niifiles(x).bytes > 28800000 && niifiles(x).bytes < 29000000
                if anat_count == 1
                    fprintf('Renaming %s to %s\n',niifiles(x).name,'anat.nii and compressing...');
                    movefile(fullfile(rawsid,niifiles(x).name),fullfile(rawsid,'anat.nii'));
                    gzip(fullfile(rawsid,'anat.nii'));
                    delete(fullfile(rawsid,'anat.nii'));
                    anat_count = anat_count + 1;
                elseif anat_count == 2
                    fprintf('Renaming %s to %s\n',niifiles(x).name,'anat_2.nii and compressing...');
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
                case strcmp(niifile{1},sprintf('run%d.nii',epi_count))
                    newname = sprintf('epi_r%.2d.nii', epi_count);
                    epi_count = epi_count + 1;
                case strcmp(niifile{1},sprintf('rest%d.nii',rest_count))
                    newname = sprintf('rest_r%.2d.nii', rest_count);
                    rest_count = rest_count + 1;
                case strcmp(niifile{1},'MPRAGE.nii')
                    newname = 'anat.nii';
                case strcmp(niifile{1},'DTI_64.nii')
                    newname = 'dti.nii';                 
                case strcmp(niifile{1},'fieldmap.nii')
                    newname = 'epi_fieldmap.nii';
                case strcmp(niifile{1},'fieldmap_phase.nii')
                    newname = 'epi_fieldmap_phase.nii';        
                otherwise
                    % cleanup tmp dir
                    rmdir(subjpath,'s');
                    error('unknown filetype: %s. Aborting...', niifile{1})
            end
            fprintf('Renaming %s to %s and compressing...\n',niifile{1},newname);
            movefile(fullfile(rawsid,niifile{1}),fullfile(rawsid,newname));
            gzip(fullfile(rawsid,newname));
            delete(fullfile(rawsid,newname));                            
        end
    end
    % Delete the extracted data
    rmdir(subjpath,'s');
end

% Print completed
fprintf('Conversion to nifti on %d subjects completed...\n', length(scannerlist));
cd(root)