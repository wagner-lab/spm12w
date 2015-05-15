function niidata = spm12w_readnii(varargin)
% spm12w_readnii(niifile, range, voxels, hdr_only, vox_only)
%
% Input
% -----
% niifile  : Full path to nifti file(s) for importing. Can accept .gz files.
%            If multiple files are provided, they must be in a cell array.         
%
% range    : Range of volumes to import (e.g., [5:100]).  
%
% voxels   : Matrix of voxels in (e.g.,[x,y,z; x,y,z]) to extract data from.
%
% sphere   : Center coordinate and diameter (d) of a sphere from which to 
%            extract values (e.g., [x,y,z,d]). Data for all voxels within the
%            sphere will be returned unless vxavg = 1 (see below). 
%
% mask     : Mask image of zeros and ones from which voxel data will be
%            extracted. Must be in same space as input nifti files. Data for 
%            all voxels within the mask will be returned unless vxavg = 1.  
%
% vxavg    : Return only voxelwise average of extracted values. This is for 
%            both voxels and image masks. In other words, thsi returns the 
%            average of all voxels in an ROI defined either by voxels or a mask.
%            (default=0)
%
% hdr_only : Flag to return only the hdr and not load any data (default=0) 
%
% vox_only : Flag to return only the voxel, sphere or mask values and not the 
%            full data matrix (default=1). Only functions if voxels, spehres
%            or masks have been specified. 
%
% resample : Type of resampling for extracting voxel values. Resampling can
%            occur if voxels are off the voxel grid (i.e., if voxel step size 
%            is 3 but voxels are based on coordiantes from a different
%            space). 0 = nearest neighbor, 1 = trilinear interpolation.
%            (default = 0). 
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
% matrix and/or the data at individual voxels or within image masks. 
%
% Examples:
%
% Read in nifti data for a range of volumes
%   >> niidata = spm12w_readnii('niifiles', 'epi_r01.nii', 'range', [1:50])
%
% Extract data for a single voxel 
%   >> vxdata = spm12w_readnii('niifiles', 'epi_r01.nii', 'voxels', [10,10,-22])
%
% Extract data for multiple voxels and average the result (voxelwise)
%   >> vxdata = spm12w_readnii('niifiles', 'epi_r01.nii', ...
%                              'voxels', [10,10,-22;13,13,26;30,50,30],...
%                              'vxavg', 1)
%
% Extract data for multiple voxels and a mask in multiple files with different
% ranges and trilinear resampling
%   >> vxdata = spm12w_readnii('niifiles', {'epi_r01.nii','epi_r02.nii'}, ...
%                              'range', {[1:50],[10:60],[30:40]}, ...
%                              'voxels', [10,10,-22;13,13,26;30,50,30], ...
%                              'mask', 'shatners_bassoon.nii',...
%                              'resample', 1}
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: November, 2014 | Updated: May, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('niifiles','', 'range','', 'voxels',[], 'mask','', ...
                'vxavg',0, 'hdr_only',0, 'vox_only',1, 'resample',0,...
                'loglevel', 1);
args = spm12w_args('nargs',2, 'defaults', args_defaults, 'arguments', varargin);

% If single niifile is given as string, convert to cell array.
if ischar(args.niifiles)
    args.niifiles = cellstr(args.niifiles);
end

% If single range entered as matrix, convert to cell array
if isnumeric(args.range)
    args.range = {args.range};
end

%Print debug message to see which files are being loaded
spm12w_logger('msg','Reading in niifile(s)...','level',args.loglevel)
if args.loglevel == 1
    for niifile = args.niifiles          
            spm12w_logger('msg',sprintf('[DEBUG] Reading nifti file at:%s', ...
            niifile{1}), 'level',args.loglevel)
    end
end

%Get hdr(s) and trim if needed (make sure number of ranges = number of hdrs)
hdrs = spm_vol(args.niifiles);
if ~isempty(args.range)
    if numel(args.range) ~= numel(hdrs)
        spm12w_logger('msg',sprintf(['[EXCEPTION] Number of range elements ',...
            '(%d) does not match number of input files (%d).'],...
            numel(args.range), numel(hdrs)),'level',args.loglevel);
        error(['Number of range elements (%d) does not match number of ',...
            'files (%d). Aborting...'],numel(args.range),numel(hdrs)); 
    else      
        % Trim the hdrs according to range.
        for hdr_i = 1:numel(hdrs)
            hdrs{hdr_i} = hdrs{hdr_i}(args.range{hdr_i});
            spm12w_logger('msg',sprintf(['[DEBUG] Triming nifti file. ', ...
               'Final size: %d volumes.'],length(hdrs{hdr_i})), 'level',...
               args.loglevel)
        end
    end
end

% Iterate through all the files, extracting data and voxels as needed.
niidata = struct('hdrs',{}); % Init empty structure
for hdr = hdrs
    % Get data if hdr_only = 0, or vox_only = 0 otherwise ignore vox_only
    if args.hdr_only == 0
        if (~isempty(args.voxels) || ~isempty(args.mask) || ~isempty(args.sphere)) && args.vox_only ==0
            data = spm_read_vols(hdr{1});  
        elseif isempty(args.voxels) && isempty(args.mask) && isempty(args.sphere)
            data = spm_read_vols(hdr{1});  
        end
    end

    % Get data at voxels (note voxdata is timeXvoxels)
    if ~isempty(args.voxels)
        voxels = args.voxels;
        spm12w_logger('msg',sprintf(['[DEBUG] Extracting voxel data for %d ',...
                      'voxels.'],size(voxels,1)),'level',args.loglevel)          
        % Convert MNI XYZ coordinates to voxel index (vox is coordinateXvoxel)
        vox = inv(hdr{1}(1).mat) * [voxels'; ones(1,size(voxels,1))];
        % Init voxdata
        voxdata = zeros(length(hdr{1}),size(voxels,1));
        % Extract voxel values
        for vol_i = 1:length(hdr{1})
            voxdata(vol_i,:) = spm_sample_vol(hdr{1}(vol_i),vox(1,:),vox(2,:),vox(3,:),...
                                 args.resample);
        end
        % Average voxdata if user requests (i.e., avg in roi)
        if args.vxavg == 1
            voxdata = nanmean(voxdata,2);
        end
    end

    % Get data at sphere
    if ~isempty(args.sphere)
        %Load mask, find where mask = 1 and apply to Y
        mask_hd  = spm_vol(mask);
        mask_vol = spm_read_vols(mask_hd);
        mask_fname = spm_str_manip(mask_hd.fname,'t');
        mask_vx  = length(find(mask_vol));
    use spm ROI
    
    end
    
    % Get data at masks
    if ~isempty(args.mask)
        %Load mask, find where mask = 1 and apply to Y
        mask_hd  = spm_vol(mask);
        mask_vol = spm_read_vols(mask_hd);
        mask_fname = spm_str_manip(mask_hd.fname,'t');
        mask_vx  = length(find(mask_vol));
    use spm ROI
    
    end
    

    
    
    % Create return structure
    niidata(end+1).hdrs = hdr{1};
    for datavars = {'data','voxels','voxdata','maskdata'}
        if exist(datavars{1}, 'var')
            niidata(end).(datavars{1}) = eval(datavars{1});
        end   
    end
end
