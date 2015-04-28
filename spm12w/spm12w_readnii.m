function niidata = spm12w_readnii(varargin)
% spm12w_readnii(niifile, range, voxels, hdr_only, vox_only)
%
% Input
% -----
% niifile  : Full path to nifti file for importing. Can accept .gz files.  
%
% range    : Range of volumes to import (e.g., [5:100]).  
%
% voxels   : Matrix of voxels in (e.g.,[x,y,z; x,y,z]) to extract data from.
%
% hdr_only : Flag to return only the hdr and not load any data (default=0) 
%
% vox_only : Flag to return only the voxel values and not the full data
%            matrix or header (default=0)
%
% loglevel : Optional log level for spm12w_logger corresponding (default=1)
%
% Returns
% -------
% niidata : A structure containing the fields hdr, data, and voxels. If
%           hdr_only=1 then only the hdr is returned if vox_only=1 then only 
%           the data for the extracted voxels is returned along with the voxels 
%           field which lists the coordiantes and the hdr. 
%
% All purpose tool for reading in nifti data. Returns the header, the full data
% matrix and/or the data at individual voxels. Call with property name/value 
% pairs.
%
% Examples:
%
% Read in nifti data for a range of volumes
%   >> niidata = spm12w_readnii('niifile', 'epi_r01.nii', 'range', [1:50])
%
% Extract data for a single voxel 
%   >> vxdata = spm12w_readnii('niifile', 'epi_r01.nii', 'voxels', [10,10,-22],
%                              'vox_only', 1)
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: November, 2014 | Updated: December, 2014
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('niifile','', 'range',[] , 'voxels',[], ...
                'hdr_only',0, 'vox_only',0, 'loglevel', 1);
args = spm12w_args('nargs',2, 'defaults', args_defaults, 'arguments', varargin);

% Get hdr and trim if needed 
spm12w_logger('msg',sprintf('Reading niifile: %s', ...
              spm_str_manip(args.niifile, 't')), 'level',args.loglevel)
hdr = spm_vol(args.niifile);
if ~isempty(args.range)
    hdr = hdr(args.range);
end

% Get data
if args.hdr_only == 0 && args.vox_only == 0
    data = spm_read_vols(hdr);    
end

% Get data at voxels
if ~isempty(args.voxels)
    voxels = args.voxels;
    voxdata = zeros(length(hdr),size(voxels,1));
    % Iterate over voxels
    for ivox = 1:size(voxels,1)
        spm12w_logger('msg',sprintf(['[DEBUG] Extracting voxel data at: ' ...
            '%d,%d,%d'], voxels(ivox,:)), 'level',args.loglevel)
        % Convert XYZ to vx index and extact timeseries
        vox = inv(hdr(1).mat) * [voxels(ivox,:) 1]';   
        % Extract timeseries 
        for ii = 1:length(voxdata)
            voxdata(ii, ivox) = spm_sample_vol(hdr(ii),vox(1),vox(2),vox(3),1);
        end      
    end
end

% Create return structure
niidata.hdr = hdr;
for datavars = {'data','voxels', 'voxdata'}
    if exist(datavars{1}, 'var')
        niidata.(datavars{1}) = eval(datavars{1});
    end   
end