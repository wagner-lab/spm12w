function spm12w_prepare(varargin)
% spm12w_prepare(rawformat)
%
% Inputs
% ------
% rawformat:    Format of raw data to prepare. Options are 'nifti', 
%               'parrec' (Philips) or 'dicom' (Siemens). Default: 'parrec'
%
% Converts either raw philips format par/rec files into nifti format files
% or raw nifti files archived from the scanner, renaming, anonymizing and
% copying the data into the study's raw directory. 
%
% spm12w_prepare should be run from the study's root directory, where it
% will look for a dbic_to_subid.txt file that specifies the correspondance
% between DBIC IDs (i.e., the name of the archive of raw data) and subject
% IDs (i.e., the anonymized name of the subject for the current study). In
% additions spm12w_prepare expects the raw data to be stored in
% raw_dbicID.tar.gz files for par/rec or nifti_dbicID.tar.gz files for
% scanner nifti files. These can be made readily on linux or you can use 
% spm12w_rawtar.m to create these files from copies of the raw data.
%
% spm12w_prepare will prepare every file it finds in the study's arch/nifti
% or arch/parrec (depending on the requested format) that is also listed
% in the dbic_to_subid.txt file. If a subject has already been converted
% (i.e., exists in the raw directory) then it will be overwritten. 
%
% spm12w_prepare workflow:
%   - Unarchive the raw data (par/rec or nifti)
%   - Convert to nifti (if par/rec)
%   - Copy to raw directory using the subject ID for the directory name.
%   - Detect anatomy and epi files and rename them from DBIC conventions to
%     spm12w conventions.
%
% Additional details on par/rec conversion: 
% PAR/REC data from 2011 stored data differently than 2012 namely, order of 
% slices volumes in the header are different. Thus neither dcm2nii nor the 
% parrec2nii command line tool work. R2AGUI appears to do the job. Also note
% that the parrec2nii par/rec converter gets the anat (or epi) transformation 
% wrong so that they are in poor alignment. r2agui seems to work better. 
% -ddw 2012
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: March, 2013 | Updated: April, 2016
% =======1=========2=========3=========4=========5=========6=========7=========8

% Input check
switch (nargin)
    case 0 
        rawformat = 'parrec';
    case 1
        rawformat = varargin{1};
    otherwise
        error('You''ve specified too many inputs. Please list only rawformat');
end

% Paths
root = pwd;
archdir = fullfile(root, 'arch');
sidfile = fullfile(archdir, 'subid_list.txt');
rawdir = fullfile(root, 'raw');

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

% Assign datadir based on rawformat
if strcmpi(rawformat, 'parrec')
    datadir = fullfile(archdir, 'parrec');
elseif strcmpi(rawformat, 'nifti')
    datadir = fullfile(archdir, 'nifti');  
else
    error('Unrecognized rawformat: %s', rawformat)
end

% Check that we are in the study's root directory and that files exists
if ~exist(sidfile, 'file')
    error(['I can''t seem to find the dbic_to_subid.txt file in the arch '...
        'directory. Are you sure it exists and that you are in the study ' ...
        'root directory? (i.e.: %s)'], root)
end

% Load subjects and sids from dbic_to_subid.txt file
filename = sidfile;
formatSpec = '%s%s%[^\n\r]';
fid = fopen(filename, 'r');
subidArray = textscan(fid, formatSpec, 'Delimiter', ',');
fclose(fid);

% Set subjects and sids to approriate entries in subidArray
subjects = deblank(subidArray{1});
sids = deblank(subidArray{2});

% Prepare either parrec or nifti files
for i = 1:length(subjects)
    cd(datadir)
    subjpath = fullfile(datadir,subjects{i});
    if strcmpi(rawformat, 'nifti') 
        % Extract the nifti files
        fprintf(['Subject: %s (%s) | Extracting nifti data from archive, ', ...
                'this may take awhile...\n'], subjects{i}, sids{i});
        untar(['nifti_',subjects{i},'.tar.gz']);
        % Find the raw NIFTI files for this subject
        files = dir(fullfile(subjpath,'*.nii'));
    else        
        % Extract the raw files
        fprintf(['Subject: %s (%s) | Extracting raw data from archive, ', ...
                'this may take awhile...\n'], subjects{i}, sids{i});
        untar(['raw_',subjects{i},'.tar.gz']);
        % Find the raw PAR files for this subject
        files = dir(fullfile(subjpath,'*.PAR'));
        options.pathpar = [subjpath,filesep];
    end
    count = 1;
    for ii = 1:length(files)
        % Prune scout scans
        if files(ii).bytes > 10000 %scouts tend to be 8078 bytes
            filelist{count} = files(ii).name;
            count = count + 1;
        end
    end
    
    % Convert from parrec if rawformat is parrec
    if strcmpi(rawformat,'parrec') 
        fprintf('Converting...\n');
        convert_r2a(filelist,options);
        fprintf('Done.\n');
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
        if strcmpi(rawformat, 'parrec')     
            moveme = fullfile(datadir,subjects{i},file_noext);
            movefile([moveme,filesep,'*.nii'],rawsid);       
        else
            moveme = fullfile(datadir,subjects{i},[file_noext,'.nii']);
            movefile(moveme,rawsid);           
        end
    end
    
    % Rename files
    if strcmpi(rawformat, 'parrec') 
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
            if exist(fullfile([subjpath,'_anat'],'anat.nii.gz'),'file')
                fprintf('Copying anonymized anatomical...\n')
                copyfile(fullfile([subjpath,'_anat'],'anat.nii.gz'), ...
                    fullfile(rawsid))
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
    else
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
    end
    
    % Delete the extracted data
    rmdir(subjpath,'s');
    clear filelist newfiles files
end

% Print completed
fprintf('Conversion to nifti on %d subjects completed...\n', length(subjects));
cd(root)