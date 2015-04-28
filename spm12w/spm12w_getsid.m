function subjects = spm12w_getsid(varargin)
% spm12w_getsid()
%
% Returns
% -------
% subjects: Cell array of subject ids.
% 
% Fetches subjects using spm_select ui. First looks in prep, if there are 
% multiple preprocessing directories then user can select which to work
% from. If there is only one preprocessing directory, then spm12w_getsid
% will immediately search in that directory. If there are no preprocessing 
% directories, spm12w_getsid will search in the raw directory. If a directory
% path is supplied. spm12w_getsid will search that location.
%
% Example:
%
%   >> spm12w_getsid()
%   >> spm12w_getsid('prep/standard')
%   >> spm12w_getsid(glm.prepdir) % If the glm parameters are in workspace
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: September 2012 | Updated: April, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Input checks. 
%       Here we do a little logic to make life easier for the user.
%       If there is a 'prep' dir, we check to see if there is more than one
%       preprocessing dir. 
%       If there is not, then we start in the only preprocessing dir. 
%       If there is, then we start in 'prep' and let the user decide
%       which preprocessing dir to work with. 
%       If there is no preprocessing dirs, we get sids from the raw directory.
switch nargin  
    case 0 
        if exist(fullfile(pwd,'prep'),'dir')
            [~,subdirs] = spm_select('List',fullfile(pwd,'prep'));
            subdirs = cellstr(subdirs);
            if length(subdirs) == 1
                startdir = fullfile(pwd,'prep',subdirs{1});
            else
                startdir = fullfile(pwd,'prep');
            end
        elseif exist(fullfile(pwd,'raw'),'dir')
            startdir = fullfile(pwd,'raw');
        else
            error(['There is no viable starting directory (i.e., prep or raw)'...
                ' from which to look for subjects. Aborting...'])
        end       
    case 1
        startdir = varargin{1};      
    otherwise
        error(['Too many arguments, spm12w_getsub takes only a starting '...
            'directory'])
end

% Select subjects starting in startdir
subjects_tmp = spm_select(Inf,'dir','Please select subjects(s)...',[],startdir);
subjects_tmp = spm_file(subjects_tmp, 'basename');  % Remove path and keep sid

% Check for negligent user
if isempty(subjects_tmp)
    error('No subjects selected. Aborting...');
end

% Convert char array to cell array
subjects = cellstr(subjects_tmp)';