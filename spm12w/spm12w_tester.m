function spm12w_tester(varargin)
% spm12w_tester(tests, sids)
%
% Inputs
% ------
% tests: Number (e.g., 3) or vector of numbers (e.g. [1,2,3,4,5]) or the
%        special word 'all' to perform all tests. (Default: 'all')
%
% sids:  Cell array of subject ids. If left unspecified, a dialog box 
%        will appear asking the user to select subjects.
%
% parloop: Optional input (0 or 1) to enable parallel processing. 
%          Default is no parallel processing (default = 0). 
%
% parcores: Optional input to specify the size of the parallel processing pool.
%           Default value is 4 cpu cores but can range from 2 to 64. Caution
%           should be used when enabling this on a shared server environment.
%
% This testing function is designed to run a new installation of spm12w
% through most of its core functions to ensure everything is functioning 
% properly. It is also useful for ensuring that new changes to the codebase 
% do not break the pipeline. This function will test the example dataset 
% (i.e., tutorial_data) on all of the following tests:
%
% 0  = Prepare the data (i.e., extract archived data and copy & rename to raw)
% 1  = Preprocess the data (normalize to epi template)
% 2  = Preprocess the data (spm12 segmentation based normalization)
% 3  = Preprocess the data (spm12 dartel normalization)
% 4  = 1st level glm build and compute.
% 5  = 1st level contrasts of task effects
% 6  = 2nd level random effects analysis
% 7  = ROI analysis
% 8  = VOI timeseries extraction
% 9  = PPI regressor creation
% 10 = Reset tutorial data to original state (warning: will viciously delete files)
%
% spm12w_tester must be run from the example dataset root directory...
%
% Examples:
%
% To run spm12w_tester on tests 1-3 on subjects s01 and s02:
%  
%       >> spm12w_tester('tests',[1,2,3],'sids',{'s01','s02'})
%
% To run spm12w_tester on tests 1,4,& 5 with parallel processing
%  
%       >> spm12w_tester('tests',[1,4:5],'sids',{'s01','s02','s03'}, ...
%          'parloop',1,'parcores',3)
%
% To run spm12w_tester on all tests while manually selecting subjects, call
% it without argument:
%    
%       >> spm12w_tester
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: February, 2013 | Updated: March, 2016
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('tests','all', 'sids','','parloop',0, 'parcores', 4);
args = spm12w_args('nargs',0, 'defaults', args_defaults, 'arguments', varargin);

% If tests is keyword 'all' set tests to tests 1 to 9
if strcmpi(args.tests,'all')
    args.tests = 1:9;
end

% If tests contains test 0 (prepare the data) do that first if sids is
% empty otherwise there are no sids to select for other tests.
if any(args.tests==0)
    % Prepare 
    % Extract, copy and rename files according to spm12w convenctions.
    % Copied data will be in root/raw directory.
    spm12w_logger('msg',sprintf(['Data is not prepared... Preparing ', ...
          'subject selection']))    
    spm12w_prepare
end

% If sids is empty and user is requesting preprocessing, check in raw instead of
% prep. 
if isempty(args.sids)     
    if any(ismember(args.tests, [1,2,3])) %For test 1:3 check in raw
        args.sids = spm12w_getsid(fullfile(pwd,'raw'));
    elseif any(ismember(args.tests,[2:5,7,8,9]))
        args.sids = spm12w_getsid;
    end
end

% Run tests by iterating through tests in args.tests
spm12w_logger('msg',sprintf(['Beginning tutorial dataset tests on %d ', ...
              'subjects'], length(args.sids)))

for testcase = args.tests
    switch testcase     
        case 1 
            % Preprocessing
            % We only need to preprocess functional datasets once for GLM &
            % Dartel. Dartel will pick up after basic prepro is complete.
            spm12w('stage','prep','sids',args.sids, 'para_file',...
                   'p_tutorial.m','parloop',args.parloop,...
                   'parcores',args.parcores);
               
        case 2 
            % Segmentation of anatomical images using p_tutorial_dartel.m
            for sid = args.sids
                spm12w_seg8(sid{1},'p_tutorial_dartel.m');     
            end
     
        case 3 % DARTEL
            % Dartel template creation followed by Dartel normalization.
            spm12w_dartel('template',args.sids,fullfile('scripts',...
                          'p_tutorial_dartel.m'));
            for i = 1:length(subjects) %TODO Double check this, it should noot loop over subjects? Or should it?
                % Anatomical to MNI
                spm12w_dartel('mprage2mni',subjects,fullfile('scripts',...
                              'p_tutorial_dartel.m'));
                % EPI to MNI
                spm12w_dartel('epi2mni',subjects,fullfile('scripts',...
                              'p_tutorial_dartel.m'));
            end
            
        case 4
            % 1st level GLM 
            % GLM will be built based on parameters file then estimated
            spm12w('stage','glm','sids',args.sids, 'para_file',...
                   'glm_tutorial.m','parloop',args.parloop,...
                   'parcores',args.parcores);

        case 5
            % 1st level contrasts
            % Using the GLM estimated in the previous step, a series of
            % contrast vectors will generated from the parameters file
            % and computed. 
            spm12w('stage','con','sids',args.sids, 'para_file',...
                   'glm_tutorial.m','parloop',args.parloop,...
                   'parcores',args.parcores);

        case 6 
            % 2nd level rfx analysis 
            % Based on the parameters file a number of 1st level contrasts
            % are brought to the 2nd level and a one-sample t-test
            % performed across conditions.
            spm12w_glm_rfx('glm_file','glm_tutorial.m','sids',args.sids);
                      
        case 7 
            % Regions of interest analysis
            % Using a combination of spherical and anatomical ROIs
            % mean parameter estimates are extracted from ROIs and submitted 
            % to basic analyses (t-test, correltions, descriptives).
            spm12w_roitool('roi_file','roi_tutorial.m','sids',args.sids);
        
        case 8 
            % VOI timeseries extraction
            spm12w_voitool('voi', subjects, fullfile('scripts','voi_tutorial.m'));
        
        case 9 
            % PPI regressor creation
            spm12w_voitool('ppi', subjects, fullfile('scripts','voi_tutorial.m'));
            
        case 10 
            % Reset the universe
            % Delete dirs and reset the tutorial dataset to a virgin state
            fprintf('Reseting example data back to original state...\n');
            % Turn off diary as it will prevent dir removal if it was left on
            diary off
            deldirs = {'qa','analysis','aux/onsets/ppi','notes','prep','raw'};
            for deldir = deldirs
                if exist(fullfile(pwd,deldir{1}),'dir')
                    fprintf('Removing directory: %s\n', deldir{1});
                    rmdir(fullfile(pwd,deldir{1}),'s')
                end
            end
            % Remove root files (except for readme.md)
            delfiles = dir(fullfile(pwd));
            for del_i = 1:length(delfiles)
                if ~delfiles(del_i).isdir && ~strcmp(delfiles(del_i).name,'readme.md')
                    fprintf('Removing file: %s\n', delfiles(del_i).name);
                    delete(fullfile(pwd,delfiles(del_i).name))
                end           
            end
    end
end

spm12w_logger('msg','Tutorial tests complete...')