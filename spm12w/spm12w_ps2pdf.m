function spm12w_ps2pdf(varargin)
% spm12w_ps2pdf('sub_id','glm_file')
%
% Inputs
% ------
% ps_file:  Path to a postscript file to convert to pdf
%
% pdf_file: Name of pdf file 
%
% spm12w_ps2pdf will attempt to determine the location of ghostscript and
% then provide this to ps2pdf.m in order to conver the ps file to pdf.
%
% Starting with Matlab 2015b, Matlab no longer ships a ghostscript mex, thereby
% breaking all our toys and causing children to cry. To ger around this we
% ask users to install GS locally and provide us with the path (by editing
% spm12path.m). Now the users are crying as well... Thanks a lot Mathworks!
%
% Examples:
%
%       >>spm12w_ps2pdf('ps_file','./prep/standard_prep/s01/preprocess.ps')
%       >>spm12w_ps2pdf('ps_file','./prep/standard_prep/s01/preprocess.ps',...
%                       'pdf_file','./output.pdf')
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: December, 2015 | Updated: December, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('ps_file','', 'pdf_file','');
args = spm12w_args('nargs',2, 'defaults', args_defaults, 'arguments', varargin);

% First let ps2pdf search for the Matlab ghostscript MEX. If that fails,
% check for a local installation of GS in the path, if one cannot be found,
% look in some common locations. If it still cannot be found, inform the
% user that the ps file could not be converted.
try
    % This will use attempt to use the default matlab supplied GS. 
    ps2pdf('psfile',args.ps_file, 'pdffile',args.pdf_file, ...
           'gspapersize','a4','deletepsfile',1);
catch
    % If on Matlab 2015b or higher, we need to find a local installation of
    % ghostscript. Ideally supplied by the user but if not we'll try and
    % find it ourselves.
    % Check for ghostcript on path (pc, mac, linux)
    gs_exec = {'gswin32c.exe','gswin64c.exe','gs'};
    for gsexec = gs_exec
        gspath = which(gsexec{1});
        if ~isempty(gspath)
            break
        end
    end  
    % If gs was not on path, search for it.
    if isempty(gspath)
        % Search for gspath in common places
        if ispc
            gs_loc = {'C:\Program Files\gs\','C:\Program Files (x86)\gs\'};
            gs_exec = {'gswin32c.exe','gswin64c.exe'};
            for gsloc = gs_loc
                gsdir = dir(gsloc{1});
                if ~isempty(gsdir)
                    for gsexec = gs_exec
                        if exist(fullfile(gsloc{1},gsdir(3).name,'bin',gsexec{1}),'file')
                            gspath = fullfile(gsloc{1},gsdir(3).name,'bin',gsexec{1});
                            break
                        end
                    end
                end
            end  
        else
            gs_loc = {'/usr/bin/gs','/usr/local/bin/gs'};
            for gsloc = gs_loc
                if exist(fullfile(gsloc{1}))
                    gspath = fullfile(gsloc{1});
                end
            end           
        end
    end
    
    % if gs path was found by the above steps, run ps2pdf with gscommand switch
    if ~isempty(gspath)
        % Intuit lib paths from gs executable path
        gslibpath = fullfile(gspath(1:strfind(gspath,'bin')-1),'lib');      
        ps2pdf('psfile',args.ps_file, 'pdffile',args.pdf_file, ...
               'gspapersize','a4','deletepsfile',1,...
               'gscommand', gspath, ...
               'gsfontpath', gslibpath , ...
               'gslibpath', gslibpath);  
    else    
        spm12w_logger('msg',sprintf(['[WARNING] Unable to conver ps file ',...
                  'preprocess.ps to preprocess.pdf. \nIf you are using ',... 
                  'Matlab 2015b or later you must specify a path to a ',... 
                  'local ghostscript installation in the spm12path.m file']))          
    end 
end
