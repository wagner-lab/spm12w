function spm12w_glm_contrast(varargin)
% spm12w_glm_contrast('sub_id','glm_file')
%
% Inputs
% ------
% sid:      Subject ID of subject for glm computation (e.g., 's01')
%
% glm_file: File specifying the parameters for glm design and contrasts
%           (e.g., 'glm_tutorial.m'). If the path is left unspecified,
%           spm12w_glm will look in the scripts directory. <optional>
%
% spm12w_glm_contrasts calculates contrasts of weighted parameter estimates
% as defined in the con structure of the supplied glm parameters file.
%
% spm12w_glm_contrasts has a number of built in functions to make
% specifying contrasts easier. Contrasts can be specified either numerically 
% (e.g., 1 -1 0 0) or by referring to the condition names (e.g., 
% human - animal). Any unspecified condition names will be given a zero weight.
% In addition, the following keywords will automatically generate a number of 
% predefined contrasts:
%
%   'housewine' : Creates contrasts of each condition vs. baseline
%                 along with a contrast of all conditions vs. baseline.
%
%   'fir_bins'  : For FIR designs, creates contrasts of each FIR time bin
%                 for each condition vs. baseline.
%
%   'fir_hrf'   : For FIR designs, creates an HRF weighted contrast across
%                 each conditions timebin.
%
% The first argument is a sid. The second argument is the name of a glm 
% parameter files (e.g., glm_tutorial.m) and is optional. If the parameter 
% file is unspecified, matlab will prompt the user to select the file.
%
% Contrast results will be written to the glm directory specified in the
% glm parameters file. 
%
% Examples:
%
%       >>spm12w_glm_compute('sid','s01')
%       >>spm12w_glm_compute('sid','s01','glm_file', ...
%                           './scripts/username/glm_tutorial.m')
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: February, 2015 | Updated: February, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('sid','', 'glm_file','');
args = spm12w_args('nargs',2, 'defaults', args_defaults, 'arguments', varargin);

% Load glm parameters
glm_con = spm12w_getp('type','glm', 'sid',args.sid, 'para_file',args.glm_file);

% Setup directories for Contrast
spm12w_dirsetup('dirtype','con','params',glm_con);

% Goto glm directory
cd(glm_con.glmdir)

% Load glm structure from estimation (could also load the SPM.mat)
spm12w_logger('msg',sprintf(['Loading contrasts from: ', ...
              '%s.mat'],glm_con.glm_name),'level',glm_con.loglevel)
load(sprintf('%s.mat',glm_con.glm_name));

% Assign con structure to glm variable.
% This is in case the con structure has changed since estimation (i.e.,
% user edited the contrasts in the glm parameters file after estimation). 
glm.con = glm_con.con;

% Reset xCon to only keep the effects of interest contrast
if length(glm.SPM.xCon) > 1
    spm12w_logger('msg','[DEBUG] Resetting previous xCon structure...', ...
                  'level',glm_con.loglevel)
end
glm.SPM.xCon=glm.SPM.xCon(1);   

% Parse the GLM model to get parameters.
% There's a number of ways to do this, at present I think the best way is
% to use the actual regressor names in SPM.xX and then match those rather
% than parse through SPM.Sess.U and SPM.Sess.C and try and figure out
% events, parametrics and user regressors and their placement in the
% actual design. Unfortunately, SPM sullies the names (adds prefixes and
% suffixes) so it's not as clean... 
reg_total = size(glm.SPM.xX.X,2);
reg_nuiss = sum(cell2mat(strfind(glm.SPM.Sess.C.name,'r-')))+1;
spm12w_logger('msg',sprintf(['GLM Design size: %d | Number of nuissance ', ...
              'predictors: %d'],reg_total,reg_nuiss),'level',glm_con.loglevel)
% Cleanup reg_names so we can do exact matches... This might need work if
% people start adding different basis functions
reg_names = glm.SPM.xX.name;
delstrings = {'Sn(1) ','*bf(1)','^1'}; %user better never use these
for delme = delstrings
   reg_names = strrep(reg_names, delme{1},'');
end

% Generate a predictor number name string to show user the design
namestr = '';
for n_i = 1:length(reg_names)
    namestr = [namestr, sprintf('(%d)%s, ',n_i,reg_names{n_i})];
end
spm12w_logger('msg',sprintf('GLM Design: %s',namestr(1:end-2)),...
              'level',glm_con.loglevel)
          
% Collapse all user conditions into a single cell array
u_conds = [glm.events, glm.blocks, glm.regressors];

% This next bit is a total hack. Because spm changes underscores to hyphens
% we need to tweak the u_conds to match so that the searches below
% proceed without error.
u_conds = strrep(u_conds,'_','-');

% Check if user requested the house wine.
% NB: The house wine always goes first in SPM.xCon. 
if isfield(glm.con,'housewine')
    spm12w_logger('msg','Pouring a glass of house wine...', ...
                  'level',glm_con.loglevel)
    % Pour user a glass of allVSbaseline (skip pmods)
    contrasts.name{1} = 'allVSbaseline';
    contrasts.vals{1} = +ismember(reg_names,u_conds); %+ to switch off logical
    spm12w_logger('msg',sprintf(['[DEBUG] Generating housewine ',...
                         'contrast(%d) %s'],length(contrasts.name), ...
                         contrasts.name{end}), 'level',glm_con.loglevel)
    spm12w_logger('msg',sprintf('[DEBUG] housewine contrast(%d) %s: %s',...
                         length(contrasts.name), contrasts.name{end}, ...
                         mat2str(contrasts.vals{end})), 'level',glm_con.loglevel)
    % Pour user a glass of each condition VSbaseline (including pmods)
    for cond = [u_conds, glm.parametrics]
        cond_index = find(strcmp(reg_names,cond{1}));
        for con_i = cond_index           
            contrasts.name{end+1} = sprintf('%sVSbaseline',reg_names{con_i});
            defcon = zeros(1,length(reg_names));
            defcon(con_i) = 1;
            contrasts.vals{end+1} = defcon;
            spm12w_logger('msg',sprintf(['[DEBUG] Generating housewine ',...
                         'contrast(%d) %s'],length(contrasts.name), ...
                         contrasts.name{end}), 'level',glm_con.loglevel)
            spm12w_logger('msg',sprintf('[DEBUG] housewine contrast(%d) %s: %s',...
                         length(contrasts.name), contrasts.name{end}, ...
                         mat2str(contrasts.vals{end})), 'level',glm_con.loglevel)
        end            
    end
end

% Parse the con structure and create complete numeric contrasts from user's 
% string, keyword and incomplete numeric contrasts. 
for cfield = fieldnames(glm.con)'
    switch cfield{1}
        case 'housewine'
            % Do nothing, housewine went first (see above)
        case 'fir_bins'
            % for now do nothing (work in progress)
        case 'fir_hrf'
            % for now do nothing (work in progress)
        otherwise
            contrasts.name{end+1} = cfield{1};
            spm12w_logger('msg',sprintf(['[DEBUG] Generating user defined ',...
                         'contrast(%d) %s'],length(contrasts.name), ...
                         contrasts.name{end}), 'level',glm_con.loglevel)
            % Parse user contrasts
            if isnumeric(glm.con.(cfield{1}))
                % If numeric, pad with zeros for lenght of design.
                lencon = length(glm.con.(cfield{1}));
                contrasts.vals{end+1} = [glm.con.(cfield{1}), ...
                                         zeros(1,reg_total-lencon)];
                % Show user the contrast values.
                spm12w_logger('msg',sprintf('[DEBUG] user defined contrast(%d) %s: %s',...
                         length(contrasts.name), contrasts.name{end}, ...
                         mat2str(contrasts.vals{end})), 'level',glm_con.loglevel)
            elseif ischar(glm.con.(cfield{1}))
                % If char, split on 'vs.' and convert to numeric
                tmpstr = strsplit(glm.con.(cfield{1}),'vs.');
                poswt = strtrim(tmpstr{1}); %remove leading/trailing whitespace
                if length(tmpstr) > 1
                    negwt = strtrim(tmpstr{2}); %remove leading/trailing whitespace
                else
                    negwt = ''; % Set to empty (i.e., no contrast of conditions)
                end
                % Now split each on wtspace to get members
                poswt = strsplit(poswt);
                negwt = strsplit(negwt);
                % Make sure to account for any underscores switched to hyphens
                poswt = strrep(poswt,'_','-');
                negwt = strrep(negwt,'_','-');
                % Convert to numeric and adjust for numel
                poswt = ismember(reg_names,poswt)/numel(poswt);
                negwt = ismember(reg_names,negwt)/numel(negwt);
                % Assign contrast vals as sum of pos and negwt
                contrasts.vals{end+1} = poswt - negwt;
                % Show user the contrast values.
                spm12w_logger('msg',sprintf('[DEBUG] user defined contrast(%d) %s: %s',...
                         length(contrasts.name), contrasts.name{end}, ...
                         mat2str(contrasts.vals{end})), 'level',glm_con.loglevel)
            else
                error('unknown data type for %s', cfield{1});
            end
    end
end

% Add contrasts to glm structure
glm.contrasts = contrasts;

% Add contrasts to prior spm structure
for con_i = 1:length(glm.contrasts.vals)
    glm.SPM.xCon(end + 1) = spm_FcUtil('Set', glm.contrasts.name{con_i}, 'T', 'c', ...
    (glm.contrasts.vals{con_i})', glm.SPM.xX.xKXs);
end

% Evaluate contrasts
spm12w_logger('msg',sprintf('Evaluating contrasts for: %s', glm_con.glm_name), ...
              'level',glm_con.loglevel)
spm_contrasts(glm.SPM);

% Save parameter structure to mat file
spm12w_logger('msg',sprintf('Saving contrast vectors to parameters file: %s.mat', ... 
               glm_con.glm_name), 'level',glm_con.loglevel)
save([glm.glm_name,'.mat'],'glm');

% Set final words (no figure, design only)
msglist{1} = glm.niceline;
msglist{2} = sprintf('GLM contrasts complete on subject: %s',glm.sid);
msglist{3} = sprintf('Parameters : %s', fullfile(glm.glmdir,[glm_con.glm_name,'.mat']));

% Print final words
for msg = msglist
    spm12w_logger('msg',msg{1},'level',glm.loglevel)
end

% return to studydir.
cd(glm.study_dir)

