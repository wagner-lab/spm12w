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
% xyzmm   : A 3 x n matrix of XYZmm coordinates (i.e., MNI). These are usually 
%           taken from the output of spm12w_readnii when using a mask or roi.
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
%   >> spm12w_writenii('niifile', 'epi_r01.nii', 'hdr', vhdr, 'data', vdata, 
%                      'xyzmm', vdata_mni)
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: November, 2014 | Updated: December, 2014
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('niifile','', 'hdr',[] , 'data',[], 'xyzmm', [], ...
                       'niizip',0, 'loglevel', 1);
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
% This is a bit inelegant for now. If no xyzmm is given then assume
% user is trying to write full data matrix and doesn't need xyzmm. 
if isempty(args.xyzmm)
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
else
    spm12w_logger('msg',sprintf('Writing niifile: %s', ...
                  spm_str_manip(args.niifile, 't')), 'level',args.loglevel)
    spm12w_logger('msg', sprintf('[DEBUG] number of vols in niifile: %d...', ... 
                  length(args.hdr)),'level',args.loglevel)
    % Figure out where voxel values go in space.
    M = args.hdr(1).mat; % voxels to mm matrix
    iM = inv(M);         % mm to voxels matrix
    DIM = args.hdr(1).dim;  %image dimensions
    [x,y,z] = ndgrid(1:DIM(1), 1:DIM(2), 1:DIM(3));    
    XYZ = [x(:)';y(:)';z(:)']; % voxel coordiantes (vx)
    clear x y z
    XYZmm = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))]; %voxel coordiantes (mm)
    clear XYZ
    % Assign input data to appropriate coordinates
    data_vol = zeros(DIM);
    for vx_i = 1:(size(args.xyzmm,2))
        idx_x = find(ismember(XYZmm(1,:),args.xyzmm(1,vx_i)));
        idx_y = find(ismember(XYZmm(2,:),args.xyzmm(2,vx_i)));
        idx_z = find(ismember(XYZmm(3,:),args.xyzmm(3,vx_i)));
        data_idx = intersect(intersect(idx_y,idx_x),idx_z);
        data_vol(data_idx) = args.data(1, vx_i);
        % Print some progress every 1000 voxels.
        if mod(vx_i,1000) == 0
            spm12w_logger('msg', sprintf('[DEBUG] %d voxels complete...', ... 
                          vx_i),'level',args.loglevel)
        end
    end
    % Not sure if we need this prep for now (stealing from my old
    % spm8w_roiwrite)
    % args.hdr = rmfield(args.hdr, 'prviate');
    % args.hdr.dt = [16 0]; %float32 precision
    % args.hdr.descrip = 'Data written by spm12w_writenii.m';
    spm_write_vol(args.hdr, data_vol);          
end



% Zip if necessary
if args.niizip 
    spm12w_logger('msg', '[DEBUG] Compressing file (gzip)...', ...
                  'level',args.loglevel)
    gzip(args.niifile);
    delete(args.niifile);
end