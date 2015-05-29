function spm12w(varargin)
% spm12w(stage, subjects)
%
% Inputs
% ------
% stage: Number (e.g., 3) or vector of numbers (e.g. [1,2,3,4,5]) for the 
%        specific analysis stage desired. Supplying a stage will override
%        override the gui and go directly to performing the analysis stage.
%
% sids: Cell array of subject ids. If left unspecified, a dialog box will
%       appear asking the user to select subjects. If stage is empty, then
%       specifying sids will not override the gui. 
%
% spm12w gui placeholder (to be finalized after everything else is working)
%
% Examples:
%
% To run spm12w preprocessing on subjects s01 and s02 overriding the gui:
%  
%       >> spm12w('stage',1,'sids',{'s01','s02'})
%
% To select analysis via the spm12w gui, type spm12w with no argument:
%    
%       >> spm12w
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: May, 2015 | Updated: May, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('stage','', 'sids','');
args = spm12w_args('nargs',0, 'defaults', args_defaults, 'arguments', varargin);

