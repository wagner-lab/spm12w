function parameters = spm12w_logger(varargin)
% spm12w_logger(msg, level, params)
%
% Input
% -----
% msg      : Message to be printed for logging. 
%
% level    : Level of log to display (default=1)
%            loglevel=1 will display all logs. 
%            loglevel=2 will display all but [DEBUG].
%            loglevel=3 will display only [INFO]
%            loglevel=4 will display no log msg whatsoever.
%
% params   : A structure of parameters (i.e., p, glm, etc.) 
%            Required only when setting up log (e,g., msg='setup_prep')
%            otherwise can be omitted.
%
% Prints a logger message with associated timestamp (in python logger
% style). If the logger message is one of the special keywords (i.e.,
% 'setup_prep' or 'setup_glm') then an associated parameters structure 
% (e.g., p or glm) must also be specified. Spm12w_logger will then 
% create the log file for the relevant state of analysis in the approrpriate
% direactory. 
%
% Setup logfile options: 
%       'setup_prep': Creates the preprocessing logfile using the name from 
%                     p.preplog.
% 
%       'setup_glm': Creates the glm logfile using the name from glm.glmlog
%
% Examples:
%
% Create the preprocessing logfile
%   >> spm12w_logger('msg', 'setup_prep', 'level', 1, 'params', p)
%
% Print a logging message 
%   >> spm12w_logger('msg', 'Beginning preprocessing on subject s01')
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: November, 2014 | Updated: December, 2014
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('msg','', 'level',1, 'params', '');
args = spm12w_args('nargs',2, 'defaults', args_defaults, 'arguments', varargin);

% Assign params to p for cleaner code. 
logp = args.params;

% Adjust internal variables for prep or glm
if strcmp(args.msg,'setup_prep')
    logpath = logp.datadir;
    logfile = logp.preplog;
elseif strcmp(args.msg,'setup_glm')
    logpath = logp.glmdir;   % this is the full path to the glmdir 
    glmname = logp.glm_name; % assign the glm name to the one the user supplied
    logfile = eval(['sprintf(''',logp.glmlog,''',glmname)']);
end

% Check msg for setup keywords
if strcmp(args.msg,'setup_prep') || strcmp(args.msg, 'setup_glm')
    % Check for directory. If it doesn't exist. Create it.
    if ~isdir(logpath)
        mkdir(logpath)
    end
    % Create the logfile/diary, if old logfile exists, delete it.
    if exist(fullfile(logpath, logfile),'file')
        fprintf('Prior logfile found and deleted...\n')
        diary off % In case the prior logfile is still open
        delete(fullfile(logpath, logfile));
    end
    fprintf('Creating logfile at: %s\n', fullfile(logpath, logfile))
    diary(fullfile(logpath, logfile))
else
    % Print message using python logger format, assume msg is INFO unless
    % using another keyword [DEBUG] or [EXCEPTION] or [WARNING].
    msg_keywords = {'[DEBUG]','[EXCEPTION]','[WARNING]'};
    if ~any(strncmpi(args.msg, msg_keywords,7)) %debug is smallest at 7 characters
        args.msg = ['[INFO] ', args.msg];   
    end
    timestamp = datestr(now, 'dd/mm/yyyy HH:MM:SS AM');
    
    switch (args.level)
        case 1
            fprintf('%s %s\n', timestamp, args.msg);
        case 2
            if ~any(strncmpi(args.msg, '[DEBUG]',7))
                fprintf('%s %s\n', timestamp, args.msg);
            end
        case 3
            if ~any(strncmpi(args.msg, msg_keywords,7))
                fprintf('%s %s\n', timestamp, args.msg);
            end
        otherwise
            % Do not print any logs for 4 or more.            
    end
end