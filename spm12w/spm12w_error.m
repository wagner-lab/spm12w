function spm12w_error(gotodir, errormsg, MError)
% spm12w_error(gotodir, errormsg, MError)
%
% Inputs
% ------
% gotodir:    Directory to cd into after exception.
%
% errormsg:   spm12w error message to display the the user
%
% MError:     MException object with information about the thrown exception
%
% spm12w_error is an internal function to alert the user to a raised
% exception and provide information to help diagnose the error.
% spm12w_error also takes an MException object which details the matlab
% specific information concerning the thrown exception. 
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: November, 2014 | Updated: December, 2014
% =======1=========2=========3=========4=========5=========6=========7=========8

% Example of how to use spm12w_error in scripts... 
% try
%   cd(p.rawdir);
% catch MError
%   errormsg = sprintf(['I can''t seem to get to the subject''s raw dir. ' ...
%                        'Are you certain your subject''s raw dir \n' ...
%                        'exists at: %s ?'], p.rawdir);
%   spm12w_error(cwd, errormsg, MError);
% end

% goto supplied directory in order to put user back to their original dir
cd(gotodir)

% display spm12w error message
fprintf('\nError!\n%s\n',errormsg);
fprintf(repmat('-',1,80));
if exist(MError,'var')
    fprintf('\n===Matlab Excepetion Information');
    fprintf('\n===Identifier: %s',MError.identifier);
    fprintf('\n===Message: %s\n',MError.message);
    fprintf(repmat('-',1,80));
end

% turn off diary in case a logfile is running
diary off