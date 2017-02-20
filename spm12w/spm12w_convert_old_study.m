function spm12w_convert_old_study(varargin)
% spm12w_readnii(old_dir, subjects, new_dir)
%
% Input
% -----
% old_dir  : Path to "old" spm99/spm2w/spm8w study directory where subject's
%            functional data are stored in old_dir/subject/FUNCTIONAL and
%            their anatomical data are stored in old_dir/subject/ANATOMY.
%
% subjects : Cell array of subject names to plunder for "old" style functional
%            (bold1_0???.img/hdr) and anatomical images (mprage.img/hdr). 
%
% new_dir  : Path to new directory where 3D (anat.nii.gz) and 4D
%            (epi_r0?.nii.gz) files will be written.
%
% Simple tool to convert 3D + time img/hdr data to 4D nifti files and 3d
% img/hdr anatomy data to 3D nifti files. This tool is designed for a specific
% use case, namely converting old study folders run with our labs' former
% spm99/spm2w/spm8w scripts to the new directory structure and file naming
% scheme of spm12w. This script takes an input directory (old_dir), a cell
% array of subject directory names and the path to a new directory to house
% the spm12w data in. Each subjects FUNCTIONAL and ANATOMY directories will be
% read and any bold1_0???.img/hdr files and mprage.img/hdr files will be
% converted to nifti and written to the new_dir under the 'raw' directory.
% In addition, old subjects will be renamed to the new sid naming scheme (i.e,
% s01,s02,s03, etc.) and a file will be writen in each directory specifying
% which of the original subjects each new sid corresponds to.
%
% This script is only really of use to users of our older spm99/spm2w/spm8w 
% scripts.
%
% Examples:
%
% Convert the following subjects un-preprocessed functional and anatomy data
% to nifty and write to a new study directory.
%   >> spm12w_convert_old_study('old_dir', '/lab/neurodata/2009_FACES', ...
%      'subjects', {'04nov09aa','04dec09bb','05dec09cc'}, ...
%      'new_dir', '/lab/neurodata/faces_study')
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: February, 2017 | Updated: February, 2017
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('old_dir','','subjects','','new_dir','');
args = spm12w_args('nargs',6, 'defaults', args_defaults, 'arguments', varargin);

% Assign output directory to raw
output_raw = fullfile(args.new_dir,'raw');

% Send message
spm12w_logger('msg',sprintf('Preparing to convert %d subjects...', ...
              length(args.subjects)));
spm12w_logger('msg',sprintf('Input directory: %s', args.old_dir));        
spm12w_logger('msg',sprintf('Output directory: %s', args.new_dir));

% Iterate over subdirs to convert and copy.
for i=1:length(args.subjects)
    spm12w_logger('msg',sprintf('Converting subject %s (%d of %d)...\n', ...
                  args.subjects{i}, i, length(args.subjects)));
    % Make sid dir
    sid_dir = fullfile(output_raw,sprintf('s%02d',i));
    if ~isdir(sid_dir)
        spm12w_logger('msg',sprintf('[DEBUG] Creating dir: %s',sid_dir));
        mkdir(sid_dir)
    end
    % Read in all runs of bold data and iterate while there are bolds to convert
    bold_exist = 1;
    bold_i = 1;
    while bold_exist == 1
        boldfiles = dir(fullfile(args.old_dir, args.subjects{i},'FUNCTIONAL', ...
                    sprintf('bold%d_*.img',bold_i)));
        if ~isempty(boldfiles)
            spm12w_logger('msg',sprintf('[DEBUG] Converting bold run %d...', ...
                          bold_i));           
            % Move into directory
            cd(fullfile(args.old_dir, args.subjects{i},'FUNCTIONAL'));           
            outfile = fullfile(sid_dir,sprintf('epi_r%02d.nii',bold_i));
            % Use spm's builtin 3d to 4d merge tool. 
            spm_file_merge({boldfiles(:).name},outfile,0);
            spm12w_logger('msg',sprintf('[DEBUG] Writing %s to:%s', ...
                          sprintf('epi_r%02d.nii.gz',bold_i),sid_dir));
            gzip(outfile)
            % gzip leaves the original behind so delete it.
            delete(outfile)
            % Also delete the junk .mat files spm_file_merge creates.
            % These aren't necessary as the .mat transforms are in the
            % nifti header.
            [tmp_dir,tmp_file] = fileparts(outfile);
            delete(fullfile(tmp_dir,sprintf('%s.mat',tmp_file))) 
            bold_i = bold_i + 1;
        else
            bold_exist = 0;
        end 
    end
    % Convert anatomy if found.
    anatdir = fullfile(args.old_dir, args.subjects{i},'ANATOMY');
    if exist(fullfile(anatdir,'mprage.img'),'file') == 2
        spm12w_logger('msg',sprintf('[DEBUG] Converting anatomy to nifti...'));
        % Convert 3d analyze to 3d nifti (can't use spm_file_merge since
        % anat is 3d). 
        cd(anatdir)
        mprageV = spm_vol(fullfile(anatdir,'mprage.img'));
        mprageI = spm_read_vols(mprageV);
        mprageV.fname = fullfile(sid_dir,'anat.nii');   
        spm_write_vol(mprageV,mprageI);
        gzip(fullfile(sid_dir,'anat.nii'));
        delete(fullfile(sid_dir,'anat.nii'));  
    end
    % Make an empty file so we now the origin of this subject and leave a trace
    % in the new study directory.
    outfile = fullfile(sid_dir,sprintf('%s_is_%s',sprintf('s%02d',i), ...
                       args.subjects{i}));
    fclose(fopen(outfile, 'w'));
    spm12w_logger('msg',sprintf('Subject: %s complete and is now %s...', ...
                  args.subjects{i}, sprintf('s%02d',i)));
end

% Final words
spm12w_logger('msg',sprintf('Finished converting %d subjects...', ...
              length(args.subjects)));

% E.T. phone home.
cd(cwd)


