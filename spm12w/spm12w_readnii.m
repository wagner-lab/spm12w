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
% sphere   : Center coordinate and radius (in mm) of a sphere from which to 
%            extract values (e.g., [x,y,z,radius]). Data for all voxels within the
%            sphere will be returned unless vox_avg = 1 (see below). 
%
% mask     : Mask image of zeros and ones from which voxel data will be
%            extracted. Must be in same space as input nifti files. Data for 
%            all voxels within the mask will be returned unless vox_avg = 1.  
%
% vox_avg    : Return only voxelwise average of extracted values. This is for 
%            both voxels and image masks. In other words, thsi returns the 
%            average of all voxels in an ROI defined either by voxels or a mask.
%            (default=0)
%
% vox_svd    : Returns first eigenvariate of extracted values. This is for 
%            both voxels, spheres and image masks. In other words, this returns
%            the average timeseries of all voxels in an ROI. Typically this 
%            should be used on timeseries data (4D files) and not con images 
%            (3D files). Note that vox_avg and vox_svd cannot both be 1. 
%            (default=0).
%
% hdr_only : Flag to return only the hdr and not load any data (default=0) 
%
% vox_only : Flag to return only the voxel, sphere or mask values and not the 
%            full data matrix (default=1). Only functions if voxels, spheres
%            or masks have been specified. 
%
% resample : Type of resampling for extracting voxel values. Resampling can
%            occur if voxels are off the voxel grid (i.e., if voxel step size 
%            is 3 but voxels are based on coordiantes from a different
%            space). 0 = nearest neighbor, 1 = trilinear interpolation.
%            (default = 0). 
%
% loglevel : Optional log level for spm12w_logger (default=1)
%
% Returns
% -------
% niidata : A structure containing fields describing the head, the voxels
%           used to extract data (for voxels, sphere and mask) as well as
%           the data extracted from those ROIs. If hdr_only=1 then only the 
%           hdr is returned if vox_only=1 then only the data for the extracted 
%           ROI is returned along with the hdr and the voxels that are 
%           contained within the ROI. 
%
% All purpose tool for reading in nifti data. Returns the header, the full data
% matrix and/or the data at voxels/spheres/masks.
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
%                              'vox_avg', 1)
%
% Extract data for multiple voxels, an 8mm sphere and a mask in multiple files 
% with different ranges and trilinear resampling
%   >> vxdata = spm12w_readnii('niifiles', {'epi_r01.nii','epi_r02.nii'}, ...
%                              'range', {[1:50],[10:60]}, ...
%                              'voxels', [10,10,-22;13,13,26;30,50,30], ...
%                              'sphere', [9,9,9,8], ...
%                              'mask', 'shatners_bassoon.nii',...
%                              'resample', 1}
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: November, 2014 | Updated: June, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('niifiles','', 'range','', 'voxels',[], 'sphere','',...
                'mask','', 'vox_avg',0, 'vox_svd',0, 'hdr_only',0, ...
                'vox_only',1, 'resample',0, 'loglevel', 1);
args = spm12w_args('nargs',2, 'defaults', args_defaults, 'arguments', varargin);

% if iscell check dimensions of niifille, if not cell then assume char.
if iscell(args.niifiles)
    [r,c] = size(args.niifiles);
    if r > c  
        % We want cell array of 1 row many columns, not vice versa.
        args.niifiles = args.niifiles';
    end
elseif ischar(args.niifiles)  % if niifile is not cell, assume str and convert
    args.niifiles = cellstr(args.niifiles);
else
    spm12w_logger('msg','[EXCEPTION] niifiles in unknown format...', ...
                  'level',args.loglevel);
    error ('Unknown format for niifiles...')
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
    
    % Assign hdr information to return structure
    niidata(end+1).hdrs = hdr{1};
    
    % Iterate over voxels, spheres and masks
    for roi_type = {'voxels','sphere','mask'}
        roi = args.(roi_type{1}); % Assign current roi_var to roi for elegance
        if ~isempty(roi)
            switch(roi_type{1})
                case 'voxels'
                    roi = roi'; % transpose since spm expects coordinateXvoxel)
                    spm12w_logger('msg',sprintf(['[DEBUG] Extracting voxel ',...
                                  'data for %d voxels.'],size(roi,2)),...
                                  'level',args.loglevel)   
             
                case 'sphere'
                    xY = struct('def','sphere','xyz',roi(1:3)',...
                                'spec',roi(4));
                    % Get voxels that correspond to that sphere in img space,
                    % need to take first index in case of 4D file. Otherwise 
                    % spm_ROI will give a structure error.
                    [~, roi, ~] = spm_ROI(xY, hdr{1}(1)); 
                    % Inform user of sphere details
                    spm12w_logger('msg',sprintf(['[DEBUG] Extracting voxel ',...
                                  'data for a sphere (%dmm radius) ',...
                                  'containing %d voxels.'], xY.spec, ...
                                  size(roi,2)),'level',args.loglevel) 
                    
                case 'mask'
                    xY = struct('def','mask','spec',roi);
                    % Get voxels that correspond to that sphere in img space,
                    % need to take first index in case of 4D file. Otherwise 
                    % spm_ROI will give a structure error.
                    [~, roi, ~] = spm_ROI(xY, hdr{1}(1));
                    % Inform user of sphere details
                    spm12w_logger('msg',sprintf(['[DEBUG] Extracting voxel ',...
                            'data for nask (%s) containing %d voxels.'], ...
                             spm_str_manip(xY.def, 't'), size(roi,2)), ...
                            'level',args.loglevel) 
            end
            % Convert XYZ coords to vox indices (voxidx is coordinateXvoxel)
            voxidx = hdr{1}(1).mat \ [roi; ones(1,size(roi,2))];
            % Init voxdata to store extracted values
            voxdata = zeros(length(hdr{1}),size(roi,2));
            % Extract voxel values
            for vol_i = 1:length(hdr{1})
                voxdata(vol_i,:) = spm_sample_vol(hdr{1}(vol_i),voxidx(1,:),...
                               voxidx(2,:),voxidx(3,:), args.resample);
            end
            % Average voxdata if user requests it (i.e., avg in roi)
            if args.vox_avg == 1
                voxdata = nanmean(voxdata,2);
                spm12w_logger('msg','[DEBUG] Averaging extracted voxel data.',...
                'level',args.loglevel) 
            end
            % Compute first eigenvariate if user requests it (i.e., svd in roi)
            if args.vox_svd == 1
                 spm12w_logger('msg','SVD DOESN''T WORK YET!!!',...
                'level',args.loglevel) 
            end
            % Assign extracted voxel data to appropriate field
            switch(roi_type{1})
                case 'voxels'
                    niidata(end).vox = roi;
                    niidata(end).voxdata = voxdata;
             
                case 'sphere'
                    niidata(end).spherevox = roi;
                    niidata(end).spheredata = voxdata;
                    
                case 'mask'
                    niidata(end).maskvox = roi;
                    niidata(end).maskdata = voxdata;
            end 
        end
    end
    % Assign the full data matrix if it exists.
    if exist('data', 'var')
        niidata(end).data = data;
    end  
end

