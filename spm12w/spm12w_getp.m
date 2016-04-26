function parameters = spm12w_getp(varargin)
% spm12w_getp(type, sid, par_file)
%
% Input
% -----
% type: The type of parameter file being loaded. Options are: p,glm,roi,voi
%       and des. Default is p. 
%
% sid: Optional input variable to specify the subject ID for the current 
%      parameters instance. This will fill in any subject specific paths
%      and return a parameters strcture that is tailored to the supplied
%      subject.
%
% para_file: Optional string describing the full or relative path to the
%            parameters file. If a filename is given without path
%            information spm12w_getp will look in the scripts directory 
%            relative to the current working directory. 
%
% Returns
% -------
% parameters: A structure whose fields are the parameters relevant to the
%             current parameters type. 
%
% Loads a parameters file of type p,d,r,glm,roi or voi and returns a
% structure of parameters along with any relevant spm12w defaults 
% (i.e., spm12w_defaults.m) and any fields that are overridden by parameters 
% in the user defaults file (spm12w_userdefaults.m) or in the parameters file
% itself (i.e., p_tutorial.m). This allows the user to override the default 
% paramters at 2 levels: a global spm12w_userdefaults.m or from within 
% the parameter files (i.e., p_tutorial.m). See spm12w_defaults.m or the user
% manual for more information on overriding default parameters.
%
% Examples:
%
% Without arguments, spm12w will prompt the user for a p file (default) and 
% return a struct of paramters.
%
%   >> p = spm12w_getp
% 
% With arguments, spm12w will return the relevant parameters file and, 
% optionally fill in the subject id for any paths that are subject specific.
%
%   >> p = spm12w_getp('type','glm', 'sid','s01', ...
%                     'para_file','scripts/glm_h8tjazz.m')
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: May, 2012 | Updated: April, 2016
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
arg_defaults = struct('type','p', 'sid','', 'para_file','');
args = spm12w_args('nargs',0, 'defaults',arg_defaults, 'arguments',varargin);

% Check for only one user directory under scripts
dirchk = dir(fullfile(pwd,'scripts'));
if length(dirchk) == 3 && dirchk(3).isdir
    para_path = fullfile(pwd,'scripts',dirchk(3).name);
else
    para_path = fullfile(pwd,'scripts');
end

if isempty(args.para_file)
    msg = 'Please select your parameters file.';
    args.para_file = spm_select(1,['^',args.type,'_.*\.m$'],msg,[], para_path);
end

% Check for indecisive user
if isempty(args.para_file)
    error('No parameters file selected. Aborting...');
end

% Load parameters file. If subj is set, fill in the sid. 
if exist(args.para_file, 'file')
    % Set subj to the appropriate structure type
    eval(sprintf('%s.(''sid'') = ''%s'';', args.type, args.sid))
    run(args.para_file);
elseif exist(fullfile(para_path,args.para_file), 'file')
    % Set subj to the appropriate structure type
    eval(sprintf('%s.(''sid'') = ''%s'';', args.type, args.sid))
    args.para_file = fullfile(para_path,args.para_file);
    run(args.para_file);
else
    error(['Unable to run supplied parameters file at (%s) '...
        'are you sure it exists?'], args.para_file)
end

% Set the internal p struct to the loaded pfile type
p_ = eval(args.type);

% Assign para file path to p struct
p_.para_file = args.para_file;

% Set potentially missing path variables to empty strings or else loading
% defaults will cause errors due to missing variables.
if strcmp(args.type,'p')
    zerofields = {'glm_name','rfx_name','ons_dir', 'roi_name'}; 
elseif strcmp(args.type,'roi') 
    zerofields = {'rfx_name','ons_dir','prep_name'};
elseif strcmp(args.type,'glm')
    zerofields = {'roi_name'}; 
elseif strcmp(args.type,'des')
    zerofields = {'rfx_name','ons_dir','prep_name','glm_name','roi_name'};
else
    zerofields = {};
end
for zfield = zerofields
    if ~isfield(p_,zfield{1})
        p_.(zfield{1}) = ''; % set field to empty string
    end
end
   
    
% Load defaults file, assign to var then delete def
def_file = which('spm12w_defaults.m');
run(def_file)     % def structure now in namespace
spm12w_def = def; % reassign def to another name
clear def         % clear def

% Load user defaults if exists, assign to var then delete def
if exist(which('spm12w_userdefaults.m'), 'file')
    def_file = which('spm12w_userdefaults.m');
    run(def_file)             % def structure now in namespace
    % check that def exists as it can be blank if userdefaults file
    % exists but is empty (i.e., userdefaults not set or commented out).
    if exist('def','var')
        spm12w_userdef = def; % reassign def to another name
        clear def             % clear def
    end
end

% Adjust parameters for defaults and alert user if defaults are overriden
% There are easier ways to do this but in order to give valid information
% to the user about which parameters are being overriden by which sources
% (e.g., parameters file or userdefautls) this more complicated solution
% was required.

% Set create field vars for comparions
def_fields = fieldnames(spm12w_def);
para_fields = fieldnames(p_);
if exist('spm12w_userdef','var')
    userdef_fields = fieldnames(spm12w_userdef);
else
    userdef_fields = {};
end

% Inform user which default variables are being overwritten and from where.
% start by iterating through all def_fields, checking for ismember of 
% parameters file and user defaults. If in user defaults, add them to p
% struture. If only in defaults, add them to p strucutre (this is how the 
% defaults get transferred from the defaults structure to the p structure). 
for i = 1:length(def_fields)
    if ismember(def_fields(i), para_fields) && ismember(def_fields(i), userdef_fields)
         fprintf(['Overriding user defaults paramater ''%s'' with the '...
                'value specified in the parameters file: %s...\n'], ...
                def_fields{i}, spm_str_manip(p_.para_file,'t'))   
    elseif ismember(def_fields(i), para_fields)
        fprintf(['Overriding spm12w defaults paramater ''%s'' with the '...
                 'value specified in the parameters file: %s...\n'], ...
                 def_fields{i}, spm_str_manip(p_.para_file,'t'))
    elseif ismember(def_fields(i), userdef_fields)
        fprintf(['Overriding spm12w defaults paramater ''%s'' with the '...
                 'value specified in the user defaults file...\n'], ...
                 def_fields{i})        
        p_.(def_fields{i}) = spm12w_userdef.(def_fields{i});
    else
        p_.(def_fields{i}) = spm12w_def.(def_fields{i});
    end
end

% Return parameters (assign internal p_ struct)
parameters = p_;
