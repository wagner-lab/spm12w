function [scannerlist, sids, rawformats, excludeseries] = spm12w_prepare_csvparser(varargin)
% spm12w_prepare_csvparser.m
%
% Inputs
% ------
% archdir: name of arch directory to search for scanner data archives
%          (default='arch').
%
% rootdir: root study directory (default=current working directory).
%
% scannerid: Filename of archived raw scanner data to be "prepared" for spm12w.
%            This is used to find the sid and rawformat for that file 
%            (optional).
%
% rawformat: Format of raw data to prepare. Options are 'nifti', 'parrec' 
%            (Philips) or 'dicom' (Siemens). If left unspecified,
%            spm12w_prepare use the format specifed in the third column of the
%            subid_list.txt (optional).    
%
% excludeseries: Optional dicom series to exclude from conversion (i.e., an
%                aborted scan due to tehnical malfunction). This requires you 
%                determine the series number to exclude using a dicom viewer
%                such as DicomBrowser. At the moment we only support excluding
%                one series per subject.
%             
%
% Returns
% -------
% scannerlist: Cell array of subject ids.
% 
% sids: sid name that the raw scanner data will be renamed to during
%       conversion. This corresponds to the subject identifier to be used
%       for this study (e.g., s01, s02, s03, etc.)
%
% rawformats: Format of raw data to prepare. Options are 'nifti', 'parrec' 
%            (Philips) or 'dicom' (Siemens). If left unspecified,
%            spm12w_prepare use the format specifed in the third column of the
%            subid_list.txt     
%
% excludeseries: The series number to exclude from conversion. 
%
% This funciton is meant to be used internally by spm12w_prepare.m and spm12w.m
% in order to parse subid_list files. It is not intended for use outside of
% these functions.
%
% Example:
%
%   >> spm12w_prepare_csvparser()
%   >> spm12w_prepare_csvparser('archdir', './arch',...
%                               'rootdir', '/data/tutorial_data',
%                               'scannerid, 'tutorial_01jan12aa_s01',...
%                               'rawformat', 'nifti')
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: June, 2016 | Updated: June, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('archdir',fullfile(pwd,'arch'),'rootdir',pwd, ...
                       'scannerid','','rawformat', '');
args = spm12w_args('nargs',0, 'defaults', args_defaults, 'arguments', varargin);

% Load the subid_list.txt or subid_list.csv or subid_list file
if exist(fullfile(args.archdir, 'subid_list.txt'),'file')
    sidfile = fullfile(args.archdir, 'subid_list.txt');
elseif exist(fullfile(args.archdir, 'subid_list.csv'),'file')
    sidfile = fullfile(args.archdir, 'subid_list');
elseif exist(fullfile(args.archdir, 'subid_list'),'file')
    sidfile = fullfile(args.archdir, 'subid_list');
else
    error(['I can''t seem to find the subid_list file in the arch '...
           'directory. Are you sure it exists and that %s is the study ' ...
           'root directory?'], args.rootdir)
end
% Load subjects and sids from subid_list.txt file
filename = sidfile;
formatSpec = '%s%s%s%[^\n\r]';
fid = fopen(filename, 'r');
subidArray = textscan(fid, formatSpec, 'Delimiter', ',');
fclose(fid);

% Set subjects, sids and rawformat to approriate entries in subidArray.
% transpose to get cell array size to match args.scannerid & args.rawformat.
scannerlist = deblank(subidArray{1})';
sids = deblank(subidArray{2})';
rawformats = deblank(subidArray{3})';
excludeseries = deblank(subidArray{4})';

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
excludeseries = excludeseries(scanneridx);
% If user supplied rawformat arg, then replace rawformat harvested from
% subid_list.txt with user supplied rawformat
if ~isempty(args.rawformat)
    rawformats(:) = {args.rawformat};
end