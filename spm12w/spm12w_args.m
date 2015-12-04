function u_args = spm12w_args(varargin)
% spm12w_args(nargs, defaults, arguments)
%
% Input
% -----
% nargs     : Minimum number of args required for nargin test (default=1).
%
% defaults  : A structure of property/value pairs of defaults
%
% arguments : user supplied argumetns (usually varargin). 
%
% Returns
% -------
% args: A structure whose fields are the properties and whose values are
%       either the defaults or the user supplied values.
%
% Helper function to parse user supplied arguments against defaults and
% return a structure of property/value pairs. Intended for internal use by
% other spm12w scripts. 
%
% Example:
%
% Get arguments
% 
%     >> args = spm12w_args('nargs',2, 'defaults', struct('niifile','', ...
%                           'subj,'s01'), 'arguments', varargin);
%
%
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: November, 2014 | Updated: November, 2014
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
if nargin < 6
    error(['Not enough arguments provided, please see help for required ' ...
           'input arguments...']);
end

% Define default args
args = struct('nargs',1, 'defaults','', 'arguments','');
options = fieldnames(args);

% Check for lack of name/value pairs
if rem(length(varargin),2) ~=0
    error('Input arguments are not property name/value pairs...')
end

% Iterate through arguments, overiding defaults.   
for argpair = reshape(varargin,2,[])
    prop_name= lower(argpair{1}); % if user supplied upper make lower 
    if any(strcmp(prop_name, options))
        args.(prop_name) = argpair{2};
    else
        error('Property %s is not a valid input property', prop_name)
    end       
end

% Parse user supplied inputs.
% Check if nargs is string (a common mistake)
if ischar(args.nargs)
    args.nargs = str2double(args.nargs);
end

% Check if sufficient arguments provided, else give error.
if length(args.arguments) < args.nargs
    error(['Not enough arguments provided, please see help for required ' ...
           'input arguments...']);
end

% Assign and parse defaults
u_args = args.defaults;
args_options = fieldnames(args.defaults);

% Check for lack of name/value pairs
if rem(length(args.arguments),2) ~=0
    error('Input arguments are not property name/value pairs...')
end

% Iterate through arguments, overiding defaults, return u_args.
for argpair = reshape(args.arguments,2,[])
    prop_name= lower(argpair{1}); % if user supplied upper make lower 
    if any(strcmp(prop_name, args_options))
        u_args.(prop_name) = argpair{2};
    else
        error('Property %s is not a valid input property', prop_name)
    end       
end