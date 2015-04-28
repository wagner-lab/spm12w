function spm12w_writenii(varargin)
% spm12w_writenii(niifile, hdr, data)
%
% Input
% -----
% niifile : Full path file name to use for writing nifti file. Will gzip if
%           '.gz' extension is found in niifile name.
%
% hdr     : Nifti header to use for writing (should be same size as data)
%
% data    : Data to write to niftifile. 
%
% loglevel : Optional log level for spm12w_logger corresponding (default=1)
%
% All purpose tool for writing nifti data. Requires a filename (can be .gz) 
% a header and data (3D or 4D). Length of header must match length of data.
% Header filenames need not match niifile variable as spm12w_writenii will
% adjust to the niifile name. 
%
% Examples:
%
% Read in nifti data for a range of volumes
%   >> spm12w_writenii('niifile', 'epi_r01.nii', 'hdr', vhdr, 'data', vdata)
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: November, 2014 | Updated: December, 2014
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('niifile','', 'hdr',[] , 'data',[], 'niizip',0, ...
                       'loglevel', 1);
args = spm12w_args('nargs', 6, 'defaults', args_defaults, 'arguments', varargin);

% Check hdr/data match
if ~length(args.hdr) == length(args.data)
    error('Size of header and size of data do not match...')
end

% Adjust niiname if gz
[fpath,fname,ext] = fileparts(args.niifile);
if strcmp(ext,'.gz')
    args.niifile = fullfile(fpath,fname);
    args.niizip = 1;
end

% Replace hdr filenames with user supplied filename
for i = 1:length(args.hdr)
    args.hdr(i).fname = args.niifile;
end

% Write data
spm12w_logger('msg',sprintf('Writing niifile: %s', ...
              spm_str_manip(args.niifile, 't')), 'level',args.loglevel)
spm12w_logger('msg', sprintf('[DEBUG] number of vols in niifile: %d...', ... 
              length(args.hdr)),'level',args.loglevel)
for i = 1:length(args.hdr)
    args.hdr(i).n = [i,1]; % adjust the volume n to reflect trimming otherwise
                           % spm_write_vol will pad the first trimmed volumes
                           % with zeros.
    spm_write_vol(args.hdr(i), args.data(:,:,:,i)); % write each vol to new file
end

% Zip if necessary
if args.niizip 
    spm12w_logger('msg', '[DEBUG] Compressing file (gzip)...', ...
                  'level',args.loglevel)
    gzip(args.niifile);
    delete(args.niifile);
end