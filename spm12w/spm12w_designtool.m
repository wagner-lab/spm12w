function spm12w_designtool(varargin) 
% spm12w_glm_build('type','','des_file','')
%
% Inputs
% ------
% type     : Type of experiment design function to run. Options are 'search',
%            'check' and 'build', (default = 'search').
%
% des_file : File specifying the parameters for design development.
%            (e.g., 'des_tutorial.m'). If the full path is left unspecified,
%             spm12w_designtool will look in the scripts directory. 
%
% spm12w_designtool performs one of three tasks depending related to designing
% fMRI experiments, checking their effeciency and finally, building preset
% trial orders out of directories of images.
% 
% The following three 'types' are available:
%
% 'search' : Design tool can SEARCH for search for an optimal design matrix
%            iterating through randomizations and null event ratios until the 
%            one with the lowest collinearity between conditions is foundes.
% 
% 'check'  : Design tool can also CHECK the collinearity on a pre-defined design
%            matrix. Useful for testing collinearity on models which require 
%            customization (multiple components to a trial, unbalanced 
%            conditions)
% 
% 'build'  : Design tool can also BUILD an fMRI design by taking as input 
%            stimulus files and outputting a numbered list of image files 
%            randomized and jittered with null event trials according to the 
%            input design matrix.
%
% Examples:
%
%       >>spm12w_designtool('type','search','des_file', ...
%                           './scripts/username/des_tutorial.m')
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: March, 2011 | Updated: April, 2016
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('type','', 'des_file','');
args = spm12w_args('nargs',2, 'defaults', args_defaults, 'arguments', varargin);

% Load glm parameters
des = spm12w_getp('type','des', 'para_file',args.des_file);

% Setup directories for GLM analysis. 
%spm12w_dirsetup('dirtype','glm','params',glm);

% Setup logfile
%spm12w_logger('msg','setup_glm', 'level',glm.loglevel, 'params',glm)

% Goto root directory
cd(des.root)

%hacks:
des.type = args.type;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check directories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make sure all directories exist unless type = 'check'
if ~strcmp(des.type,'check')
    if ~isdir(fullfile(des.root,'design'));    [s,w]=system(sprintf('mkdir "%s"',fullfile(des.root,'design')));    end
    if ~isdir(des.design_dir);                 [s,w]=system(sprintf('mkdir "%s"',des.design_dir));        end
    if ~isdir(fullfile(des.design_dir,'mat')); [s,w]=system(sprintf('mkdir "%s"',fullfile(des.design_dir,'mat'))); end
    if ~isdir(des.output_dir);                 [s,w]=system(sprintf('mkdir "%s"',des.output_dir));        end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Super switch! SEARCH, CHECK, BUILD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch(des.type)
    case{'search','SEARCH'}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Loop through des.iterate and des.nullratio  
    %%% and find best design matrix according 
    %%% to des.criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    start_time = datestr(now);   
    fprintf('\n==========Beginning design generation and selection at %s\n', start_time); 
    %%% SEARCH FOR BEST DESIGN MATRIX.
    %Check for fixed or rapid-er design
    if strcmp(des.nullratio,'fixed')
        nulls = 1;
    else
        nulls = length(des.nullratio);
    end
    for i = 1:nulls 
        effstudy = 1; effruns  = 1;
        for ii = 1:des.iterate
           %%%GENERATE MATRIX
           des.design{i}.matrix{ii} = generate_design(des.cond_specs, des.nses, des.nullbeg, des.nullend, des.nullratio, des.nullfix, des.counter);%des.nullratio(i) for normal operation april2012
           %%%REPORT STATS ON MATRIX
           actual_ratio = (sum(des.design{i}.matrix{ii}{1}==666)-des.nullbeg-des.nullend)/(length(des.design{i}.matrix{ii}{1})-des.nullbeg-des.nullend);    
           fprintf('=Block:%d Iteration:%2.0f   Run Length:%d \n=Desired Nullratio:%1.2f Actual Nullratio:%1.3f',i,ii,length(des.design{i}.matrix{ii}{1}),des.nullratio(i),actual_ratio); 
           fprintf('\n');   
           %%%GENERATE SPM DESIGN MATRIX  
           [des.design{i}.SPMstudy{ii}, des.design{i}.SPMruns{ii}] = design_maker(des.design{i}.matrix{ii}, des.cond_specs, des.nses, des.TR, des.hpf, des.autocorr);           
           %%%GENERATE AND REPORT DESIGN EFFICIENCY
           des.design{i}.effstudy{ii} = design_efficiency(des.design{i}.SPMstudy{ii});
           eff_reporter(des.design{i}.effstudy{ii}, des.design{i}.SPMstudy{ii}.Sess.U);
           if des.nses > 1
                for x = 1:length(des.design{i}.SPMruns{ii})
                    des.design{i}.effruns{ii}{x} = design_efficiency(des.design{i}.SPMruns{ii}(x));
                    maxruns(x)                 = max(des.design{i}.effruns{ii}{x});
                end
                eff_reporter(des.design{i}.effruns{ii});
           else
           fprintf('\n');   
           end
           %%%CHECK IF TOP DOG STUDY FOR STUDY AND RUNS
           if max(des.design{i}.effstudy{ii}) < effstudy
                mindex_study(i) = ii;
                effstudy        = max(des.design{i}.effstudy{ii});
           end
           if max(maxruns) < effruns
                mindex_runs(i)  = ii;
                effruns         = max(maxruns);
           end          
        end
    end
    %REPORT AND SAVE BEST DESIGNS
    for i = 1:nulls
        %%% ASSIGN WINNER ACCORDING TO des.CRITERIA
        if strcmp(des.criteria,'study')
            mindex = mindex_study(i);
        elseif strcmp(des.criteria,'runs')
            mindex = mindex_runs(i);
        else
            error('Criteria is neither study nor runs, check your spelling for des.criteria');
        end
        %GET THE RUN WIDE MIN/MAX
        for x = 1:length(des.design{i}.effruns{mindex})
            maxruns(x) = max(des.design{i}.effruns{mindex}{x});
        end
        %CALCULATE THE ACTUAL RATIO FOR THE WINNER
        actual_ratio = (sum(des.design{i}.matrix{mindex}{1}==666)-des.nullbeg-des.nullend)/(length(des.design{i}.matrix{mindex}{1})-des.nullbeg-des.nullend);
        fprintf(['\n==========Best design for nullratio of %1.2f is\n' ...
                '===Iteration #:%d Run Length:%d Actual nullratio:%1.2f\n' ...
                '===Design Wide: MinR:%1.3f | MaxR:%1.3f\n' ...
                '===Run Wide   : MinR:%1.3f | MaxR:%1.3f\n'], des.nullratio(i), mindex, ...
                length(des.design{i}.matrix{mindex}{1}), actual_ratio, ... 
                min(des.design{i}.effstudy{mindex}), max(des.design{i}.effstudy{mindex}),...
                min(maxruns), max(maxruns)); 
        %SAVE BEST DESIGN 
        dfilename = [des.design_pfx,'_',num2str(mindex),'_',num2str(des.nullratio(i)),'.mat'];
        D_RUNS    = des.design{1}.matrix{mindex};
        save (dfilename, 'D_RUNS');
        [s,w]     = system(spm8w_osbabel(sprintf('mv %s %s',dfilename,fullfile(des.design_dir,'mat'))));
        fprintf('===Design saved to file %s in directory:\n===%s\n',dfilename, [des.design_dir,'/mat']) 
        %%% Spit Time
        stop_time                 = datestr(now);
        time_elapsed              = etime(datevec(stop_time),datevec(start_time)); %time elapsed in seconds
        [hours, minutes, seconds] = spm8w_timecalc(time_elapsed);
        fprintf('==========Design generation and selection finished and took %d hours, %d minutes and %d seconds...\n',stop_time,hours,minutes,seconds);
    end

    case{'check','CHECK'}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Load a design matrix or SPM mat file 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%If SPM.mat, caculate efficiency
    %%%If design matrix, generate SPM, calculate efficiency
    dmat_file = spm_select(1,'.mat$','Please select your Design matrix or SPM file.',[],des.root);
    dmat = load(dmat_file);
    if(isfield(dmat,'SPM'))
        fprintf('\n====Design efficiency for: %s\n', spm_str_manip(dmat_file,'t')); 
        %%%GENERATE AND REPORT DESIGN EFFICIENCY
        spm_eff = design_efficiency(dmat.SPM);
        eff_reporter(spm_eff, dmat.SPM.Sess.U);
        fprintf('\n');     
    elseif(isfield(dmat,'D_RUNS'))
        %%%GATHER DEFAULTS AND GENERATE MISSING VARS 
        dmat.hpf      = 'Inf'; 
        dmat.autocorr = 'none';
        dmat.nses     = length(dmat.D_RUNS);
        dmat.TR       = spm_input('What is the TR for this study?','+1','e','2.5');
        close(gcf);
        %%%GENERATE COND_SPECS FROMS DESIGN MATRIX
        %Concat matrix
        numcond = [];
        for i=1:length(dmat.D_RUNS)
            numcond = [numcond;dmat.D_RUNS{i}];
        end
        numcond(numcond==666)=[];      %Remove nulls
        numcond = max(numcond);        %find max condition number
        %GENERATE SIMPLIFIED COND_SPECS, GIVE DUMMY NAMES TO CONDS
        for i = 1:numcond
            if i < 10
                dmat.cond_specs{i,1} = ['0',num2str(i)];
            else
                dmat.cond_specs{i,1} = num2str(i);
            end
        end
        %%%GENERATE SPM STRUCTURE
        [dmat.SPMstudy, dmat.SPMruns] = design_maker(dmat.D_RUNS, dmat.cond_specs, dmat.nses, dmat.TR, dmat.hpf, dmat.autocorr);
        fprintf('\n====Design efficiency for: %s\n', spm_str_manip(dmat_file,'t')); 
        %%%GENERATE DESIGN EFFICIENCY
        spm_eff = design_efficiency(dmat.SPMstudy);
        %%%GENERATE RUN EFFICIENCY
        for x = 1:length(dmat.SPMruns)
            spm_effruns{x} = design_efficiency(dmat.SPMruns(x));
        end
        %%%REPORT EFFICIENCY
        eff_reporter(spm_eff, dmat.SPMstudy.Sess.U);
        if length(dmat.SPMruns) > 1
            eff_reporter(spm_effruns);
        else
            fprintf('\n');
        end
    else
        error('File %s is not an appropriately formatted DM_study.mat or SPM.mat file.', spm_str_manip(dmat_file,'t')); 
    end

    case{'build','BUILD'}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Load a design matrix and rename stim
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %Load a mat file or a text file?
    dmf = spm_select(1,['^',des.design_pfx,'_.*\.mat$'],'Please select your Design matrix file.',[],des.design_dir);
    dmf = load(dmf);
    dmf = dmf.D_RUNS;
    
    %load stim filenames and add to 4th col of cond_specs      
    % Generate randomizers and add to 5th col of cond specs
    %also save onsets
    for cond_i = 1:size(des.cond_specs,1)
        tmpfiles = dir(fullfile(des.input_dir,des.cond_specs{cond_i,3}));
        des.cond_specs{cond_i,4} = {tmpfiles(~[tmpfiles.isdir]).name}';
        des.cond_specs{cond_i,5} = randperm(str2double(des.cond_specs{cond_i,2}));       
    end
    %Setup stim counter (for each condition)
    stimcounter(1:length(des.cond_specs)) = 1;
    %Rename files and output to each run
    for run_i = 1:size(dmf,2)
        fprintf('Working on run %d, size %d...\n',run_i,length(dmf{run_i}));
        % Iterate through dmf
        for trial_i = 1:length(dmf{run_i})
            if dmf{run_i}(trial_i) == 666
                infile = fullfile(des.input_dir,des.null_filename);
                condname = 'fixation';
            else
                condi = dmf{run_i}(trial_i);
                indir = fullfile(des.input_dir,des.cond_specs{condi,3});
                infile = fullfile(indir,des.cond_specs{condi,4}{des.cond_specs{condi,5}(stimcounter(condi))});
                stimcounter(condi) = stimcounter(condi)+1;     
                condname = des.cond_specs{condi,3};
            end
            outfile = fullfile(des.output_dir,sprintf('run%02d_%03d_%s.jpg',run_i, trial_i,condname));
            fprintf('Infile: %s\n',infile)
            fprintf('Outfile: %s\n',outfile)
            copyfile(infile,outfile);
        end
    end
    %save onsets
    descat = [dmf{1};dmf{2}];
    for cond_i = 1:size(des.cond_specs,1)
        con_idx = find(descat==cond_i);
        con_idx = con_idx-1;
        filename = fullfile([des.cond_specs{cond_i,1},'.txt']);
        save(filename, '-ascii', 'con_idx')        
    end
    
    
    % Previously hardcoded stuff for a specific use case. 
    
%     stimcounter(1:length(des.cond_specs)) = 1;
% 
%     for i = 1:length(dmf)  %for each run
%         randchars = randperm(8);
%         randcatch = randperm(18);
%         randnull  = randperm(6);
%         count1 = 1; count2 = 1; count3 = 1; count4 = 1; count5 = 1; count6 = 1;
%         count7 = 1; count8 = 1; count9 = 1; count10 = 1; count11 = 1; count12 = 1;
%         count13 = 1; count14 = 1; count15 = 1; count16 = 1; count17 = 1; count18 = 1;
%         count19 = 1; count20 = 1; count21 = 1; count22 = 1;
%         for ii = 1:length(dmf{i})
%             if dmf{i}(ii) == 666
%                 [s,w]=system(spm8w_osbabel(sprintf('cp "%s/"%s "%s/RUN%02d/"%04d_fixation.mp4',des.input_dir,des.null_filename,des.output_dir,i,ii)));
%             else
% %                 [namesize, tmp] = size(des.cond_specs{dmf{i}(ii),4});
% %                 if namesize > 1
% %                     if stimcounter(dmf{i}(ii)) > namesize
% %                         stimcounter(dmf{i}(ii)) = 1;
% %                         [s,w]=unix(['cp ',des.input_dir,'/',des.cond_specs{dmf{i}(ii),3}{stimcounter(dmf{i}(ii))},' ',des.output_dir,'/RUN',num2str(i),'/',num2str(ii),'_cond',num2str(dmf{i}(ii)),'.jpg']);
% %                     else
% % %                         des.cond_specs{dmf{i}(ii),3}
% % %                          des.cond_specs{dmf{i}(ii),3}{stimcounter(dmf{i}(ii))}
% %                         [s,w]=unix(['cp ',des.input_dir,'/',des.cond_specs{dmf{i}(ii),3}{stimcounter(dmf{i}(ii))},' ',des.output_dir,'/RUN',num2str(i),'/',num2str(ii),'_cond',num2str(dmf{i}(ii)),'.jpg']);
% %                     end
% %                     stimcounter(dmf{i}(ii)) = stimcounter(dmf{i}(ii))+1;
%                % else
%                    % [s,w]=unix(['cp ',des.input_dir,'/',des.cond_specs{dmf{i}(ii),4},' ',des.output_dir,'/RUN',num2str(i),'/',num2str(ii),'_cond',num2str(dmf{i}(ii)),'.jpg']);
%                    %get files in dir, randomly chose 2 and reset counter
%                    dirname = fullfile(des.input_dir,des.cond_specs{dmf{i}(ii),4});
%                    curdir  = dir([dirname,filesep,'*.mp4']);
%                    curidx  = eval(['count',num2str(dmf{i}(ii))]);
%                    if strcmp(des.cond_specs{dmf{i}(ii),4},'NULL')
%                        currand = randnull;
%                    elseif strcmp(des.cond_specs{dmf{i}(ii),4},'CATCH')
%                        currand = randcatch;
%                    else
%                        currand = randchars;
%                    end                      
%                    [s,w]=system(spm8w_osbabel(sprintf('cp "%s/"%s "%s/RUN%02d/"%04d_%s_%s',dirname,curdir(currand(curidx)).name,des.output_dir,i,ii,des.cond_specs{dmf{i}(ii),4},curdir(currand(curidx)).name)));  
%                    eval(['count',num2str(dmf{i}(ii)), ' = count',num2str(dmf{i}(ii)),' + 1;']);
%             end             
%         end
%     end
    
    
    disp('built')
    
    %Step 1-Read in design matrix. 
    %Step 2-Copy images from des.input to des.output 
    %Step 3-Output onsets.

end

%%% SOD OFF AND GO HOME
cd(des.root);


%%% ADDITIONAL FUNCTIONS
function run=generate_design(cond_specs, nses, nullbeg, nullend, nullratio, nullfix, rand_type)
% Generate a design matrix per run using paramters defined in a structure 
% usually the des. structure from D_params.m but could be used 
% seperately if neededes. Randomization is done with functions 
% taken from counter.m (unknown author). Interleaving of fixations
% is done according to nullratio. 
% MODIFICATIONS:
% -First version - DDW Mar/11
% -Included optional randomization for instances were counter.m breaks.
% -Added support for fixed designs (i.e. slow ER) - DDW April/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Init some variables
%cond = length(cond_specs);
%%% Check for daft user
%%%Starting to break some code. Commenting previous stuff - April/12 DDW
%Are all stim amounts divisible by run?
% for i = 1:size(cond_specs,1)
%     if rem(str2num(cond_specs{i,2}),nses) ~= 0
%         error('Condition stimulus number is not divisible by number of runs. Aborting....');
%     end
% end
%If using counter.m are all conditions stim number same?
if rand_type == 1
    for i = 1:size(cond_specs,1)
        tmp(i)=str2double(cond_specs{i,2});
    end
    if std(tmp) ~= 0
        error('Unbalanced stimulus numbers in conditions. Cannot use counter.m to generate design. Aborting...');
    end
clear tmp
end
%%% Pseudorandomization vs. Counterbalancings
%%% Generate matrix per run. We now use numstim/run rather than split
%%% numstim by nses. Which was problematic. There are some issues with
%%% geting the distribution even over a run but we'll figure it out as we
%%% go April/12
%%% categories are equally represented across all runs
switch(rand_type)
    case 0
    %Do randomization with constraints (pseudo "pseudo-randomization"!).
    %This can be long, but not as long as full counterbalancing.
    fprintf('\n===Generating design matrix using pseudorandomization');
    for i = 1:nses
        %%%Create design matrix for each run
        run{i} = (pseudorando(cond_specs, nses));
        fprintf('.');
    end
    
    case 1
    %Do randomization with counterbalancing
    fprintf('\n===Generating design matrix using counterbalancing');
    for i = 1:nses
        %%%Create design matrix for each run
        run{i} = (counter(length(cond_specs), (str2double(cond_specs{1,2})*length(cond_specs)/nses)))';
        fprintf('.');
    end
end
%%% ADD NULL EVENTS
%%% Interleave null events (666 The Mark of the NULL!) according to nullratio and paddes. 
if strcmp(nullratio, 'fixed')
    for i = 1:nses
        new_idx        = 1;                       %Index for new run matrix with fixations
        old_idx        = 1;                       %Index for original run matrix sans null event 
        while old_idx ~= length(run{i})+1         %iterate until we're done going through condition matrix
           %Check for duration
            dur = str2num(cond_specs{run{i}(old_idx),3});
            for ii = 1:dur
                newrun(new_idx) = run{i}(old_idx);    %use original condition number for this row.                
                new_idx = new_idx + 1;
            end
            old_idx = old_idx + 1;
            for ii = 1:nullfix
                newrun(new_idx) = 666;
                new_idx = new_idx + 1;
            end           
        end
        %Generate pading matrices
        begpad(1:nullbeg,1) = 666;
        endpad(1:nullend,1) = 666;
        run{i} = [begpad; newrun'; endpad];        %assign run matrix to the new matrix with fixations
        %CLEANUP
        clear newrun nullcontrol chancecontrol chancebroken chancecount new_idx old_idx
    end    
else
    for i = 1:nses
        new_idx        = 1;                       %Index for new run matrix with fixations
        old_idx        = 1;                       %Index for original run matrix sans null event 
        chancecontrol  = nullfix;                 %Proect against runs of 4 null events
        chancebroken   = 0;                       %she isn't broke yet.
        chancecount    = 0;
        nullcontrol    = nullratio;               %Assign nullratio to tmp var
        while old_idx ~= length(run{i})+1         %iterate until we're done going through condition matrix
            if rand() < nullcontrol               %let chance decide if this trial is a null event
                %Loop through last #chancecontrol 
                %trials to check for runs of null events
                if new_idx > 5
                    for chance_i = chancecontrol:-1:1
                        %%%If last trial and previous three were null event, chance
                        %%%is broken
                        if newrun(new_idx-chance_i) == 666 && chancecount == chancecontrol -1
                            chancebroken = 1;
                            chancecount  = 0;
                        elseif newrun(new_idx-chance_i) == 666
                        %%% If last chance_i trial was null event, add one to
                        %%% chancecount.
                            chancecount = chancecount + 1;                   
                        else
                        %%% If any one of those trials was no fixation. Chance is
                        %%% working fine.
                            chancebroken = 0;
                            chancecount = 0;
                            break
                        end
                    end
                end
                if chancebroken == 1                    %if a full run of null events was found, revert
                    newrun(new_idx) = run{i}(old_idx);  %use original condition number for this row.
                    old_idx         = old_idx + 1;      %advance old_idx for original matrix and increase
                                                        %null event ratio by 1/number of trials
                    %Not needed as we're always over anyway. -DDW mar/11
                    %nullcontrol     = nullcontrol + 1/length(run{i});   %null event odds.    
                    chancebroken    = 0;
                else 
                    newrun(new_idx) = 666;               %null event is a new condition
                end
            else
                newrun(new_idx) = run{i}(old_idx); %use original condition number for this row.
                old_idx         = old_idx + 1;     %advance old_idx
            end
            new_idx             = new_idx + 1;     %advance new matrix
        end    
        %disp(newrun')
        %Generate pading matrices
        begpad(1:nullbeg,1) = 666;
        endpad(1:nullend,1) = 666;
        run{i} = [begpad; newrun'; endpad];        %assign run matrix to the new matrix with fixations
        %CLEANUP
        clear newrun nullcontrol chancecontrol chancebroken chancecount new_idx old_idx
    end
end
%%%Drop a new line
fprintf('\n');
%%% Randomly insert null events to make all runs same length
%%% This will inflate the nullratio, but striked me as better than
%%% making all runs equal to the smallest one by randomly
%%% chopping out null events from the longer runs. 
% Determine biggest run
maxrun     = 0;
maxrun_idx = 0;
for i = 1:nses
    if length(run{i}) > maxrun
        maxrun     = length(run{i});
        maxrun_idx = i;
    end
end
for i = 1:nses
    %%% Skip the run that's already at max
    if i ~= maxrun_idx
       %%%Determine how many fixations to add
       for ii = 1:maxrun - length(run{i})
            %%%find index of place to add fixation taking into account padding
            rand_idx =  ceil(rand*(length(run{i})-nullend-nullbeg) + nullbeg); 
            run{i}   = [run{i}(1:rand_idx); 666; run{i}(rand_idx+1:end)];
       end
    end        
end
%%% CLEANUP
clear maxrun maxrun_idx i ii rand_idx
return

%%% ADDITIONAL FUNCTIONS
function [SPMstudy, SPMruns]=design_maker(matrix, cond_specs, nses, TR, hpf, autocorr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spm('Defaults','fmri');  %not sure if this is necessary right now. 
%%%---------------------------------------------------
%%% User-defined parameters for this analysis
%%%---------------------------------------------------
SPM.xY.RT          = TR;                         % experiment TR in seconds
SPM.xGX.iGXcalc    = 'None';                     % global normalization: OPTIONS:'Scaling'|'None'
SPM.xX.K.HParam    = hpf;                        % high-pass filter cutoff (secs) [Inf = no filtering] -------DDW CHANGE
SPM.xVi.form       = autocorr;                   % intrinsic autocorrelations: OPTIONS: 'none'|'AR(1) + w' ------DDW CHANGE
%%%---------------------------------------------------
%%% basis functions and timing parameters
%%%---------------------------------------------------
SPM.xBF.name       = 'hrf';         % OPTIONS:'hrf'
                                    %         'hrf (with time derivative)'
                                    %         'hrf (with time and dispersion derivatives)'
                                    %         'Fourier set'
                                    %         'Fourier set (Hanning)'
                                    %         'Gamma functions'
                                    %         'Finite Impulse Response'                          
SPM.xBF.T          = 16;            % number of time bins per scan
SPM.xBF.T0         = 1;             % reference time bin 
SPM.xBF.UNITS      = 'scans';       % OPTIONS: 'scans'|'secs' for onsets ----------------DDW CHANGE
SPM.xBF.Volterra   = 1;             % OPTIONS: 1|2 = order of convolution; 1 = no Volterra

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup design matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('==========Beginning design specification\n'); 

%%% Init Trial and User Regressors 
SPM.Sess.C.C    = [];
SPM.Sess.C.name = {};
SPM.Sess.U      = [];  

%%%PARSE EVENT-RELATED DESIGNS
%%%Runs criteria requires seperate designs per run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Concat matrix if output is from design_generator
%%% Otherwise do not concat (if input is a premade matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% if size of matrix is bigger than 1,1 then concat?
%%% Concat design matrix into one long matrix.
matrix_cat = [];
for i=1:length(matrix)
    matrix_cat = [matrix_cat;matrix{i}];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Study SPM design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Init size
SPM.nscan          = nses .* length(matrix{1});         % number of scans for each of nsess sessions
%%%Assign onsets
for i=1:size(cond_specs,1) 
    SPM.Sess.U(i).name = cond_specs(i,1);               % string in cell
    SPM.Sess.U(i).ons  = find(matrix_cat == i);         % find indices in matrix of condition number
    SPM.Sess.U(i).dur = 0;
    SPM.Sess.U(i).P.name = 'none';  
end
%%%GENERATE SPM DESIGN MATRIX
%SPM       = spm_fMRI_design(SPM);
%%%Using evalc to capture output and not polute the command window
evalc('SPM = spm_fMRI_design(SPM)');
%%%Remove the SPM.mat file that is automatically generated
%%%We can save it later.
%[s,w]     =unix('rm SPM.mat');
%%%Assign to output
SPMstudy = SPM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Runs SPM designs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Init size
for i = 1:length(matrix)
    SPM.nscan                 = length(matrix{i});            % number of scans for each run
    %%%Assign onsets
    for ii=1:size(cond_specs,1) 
        SPM.Sess.U(ii).name   = cond_specs(ii,1);             % string in cell
        SPM.Sess.U(ii).ons    = find(matrix{i} == ii);        % find indices in matrix of condition number
        SPM.Sess.U(ii).dur    = 0;
        SPM.Sess.U(ii).P.name = 'none';  
    end
    %%%GENERATE SPM DESIGN MATRIX
    %SPM       = spm_fMRI_design(SPM);
    %%%Using evalc to capture output and not polute the command window
    evalc('SPM = spm_fMRI_design(SPM)');
    %%%Remove the SPM.mat file that is automatically generated
    %%%We can save it later.
  %  [s,w]     =unix('rm SPM.mat');
    %%%Assign to output
    SPMruns(i) = SPM;
end
return

function efficiency=design_efficiency(SPM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Do pairwise correlations, saving max correlation 
%%% between predictors. Round correlations so that the max
%%% can also be negative correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Init max_r
max_r(1:length(SPM.Sess.U),1) = zeros;
for i = 1:length(SPM.Sess.U)
    for ii = 1:length(SPM.Sess.U)
        if i ~= ii
            corrtmp = abs(corrcoef(SPM.xX.X(:,i),SPM.xX.X(:,ii)));
            if max_r(i) < corrtmp
                max_r(i) = corrtmp(1,2);
            end            
        end
    end
end
efficiency = max_r;
return

function eff_reporter(eff, cond_names)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check if study-wide or run-wide efficiency and report
%%% efficiency to command window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Check input
[tmp,type] = size(eff);
if type == 1
    %%%Report Design Wide Min/Max per condition.
    fprintf('=Design wide Min/Max:\nMin:%1.3f | Max:%1.3f\n=Design wide Max per condition:\n',min(eff),max(eff)); 
    %%%Report Design Wide Min/Max per condition.
    for i = 1:length(cond_names)      %Get num conditions from SPM.Sess.U
        fprintf('%12s Max:%1.3f | ',cond_names(i).name{1},eff(i))
        if i~= 1 && i~=length(cond_names) && rem(i,3) == 0 
            fprintf('\n');  %print new line every 4 conditions
        end
    end
else
    %%%Get efficiency per run
        for i = 1:length(eff)
            maxruns(i)     = max(eff{i});
        end
    %%%Report Design Wide Min/Max per Run.
    fprintf('\n=Run wide Min/Max across conditions:\nMin:%1.3f | Max:%1.3f\n=Run wide Max across conditions:\n',min(maxruns),max(maxruns)); 
    for i = 1:length(eff)      %Get num runs
        fprintf('Run%2.0f Max:%1.3f | ',i, maxruns(i))
        if i~= 1 && rem(i,3) == 0 
            fprintf('\n');  %print new line every 4 conditions
        end
    end               
    fprintf('\n');   
end
return

function condition_matrix=pseudorando(cond_specs, nses)
%%% Randperm? 
%%% matrix(randperm(length(matrix))) = matrix
%%% Pseudo Pseudo Randomizaiton
%%% Step 1- Find greatest common factor of all stim numbers
%%%         between 1 and 4 (splits run into 1,2,3 or 4 parts).
%%%         any higher seems ridiculous.
%%% Step 2- Generate seperate matrix acording to greatest common factor.
%%% Step 3- Permute and check for two consecutive trial types. Permute
%%%         again until awesome.
%%% Step 4- Concat all matrices back together and output.
%%% FIND greatest common factor between 1 and 4.
for i = 4:-1:1
    com_factor = i;                                      %Assume i is the com_factor
    for ii = 1:size(cond_specs,1)                        %Iterate through all conditions
        if rem(str2double(cond_specs{ii,2}),i) ~= 0 %If stimnumber is not divisible by i, break
            com_factor = 0;                         
            break
        end
    end
    if com_factor == i                                   %If com_factor is found, break.
        break
    end
end
dsplit{i:com_factor} = [];
%%% GENERATE matrices of all stim, split by common factor.
for i = 1:com_factor
    dsplit{i} = [];     %init var
    for ii = 1:size(cond_specs,1)
        tmpmat(1:str2double(cond_specs{ii,2})/com_factor,1) = ii;
        dsplit{i} = [dsplit{i};  tmpmat];
        clear tmpmat
    end
end
%%% PERMUTE each matrix and check for consecutive trials. 
%%% Spacing is hardcoded to 5. (but with gonogo the spacing often has to be
%%% less than 2. 
for i = 1:com_factor
    victory = 0;
    while victory == 0
        dsplit{i} = dsplit{i}(randperm(length(dsplit{i}))); %PERMUTE
        for ii = 2:length(dsplit{i})
            %eventually figure out how to do this without harcoding using
            %des.spacing
            if ii == 2 && dsplit{i}(ii) == dsplit{i}(ii-1)
                break
%            elseif ii > 2 && dsplit{i}(ii) == dsplit{i}(ii-1) | dsplit{i}(ii) == dsplit{i}(ii-2)
%                break
%            elseif ii > 3 && dsplit{i}(ii) == dsplit{i}(ii-1) | dsplit{i}(ii) == dsplit{i}(ii-2) | dsplit{i}(ii) == dsplit{i}(ii-3)
%                break  
%            elseif ii > 4 && dsplit{i}(ii) == dsplit{i}(ii-1) | dsplit{i}(ii) == dsplit{i}(ii-2) | dsplit{i}(ii) == dsplit{i}(ii-3) | dsplit{i}(ii) == dsplit{i}(ii-4)
%                break  
%            elseif ii > 5 && dsplit{i}(ii) == dsplit{i}(ii-1) | dsplit{i}(ii) == dsplit{i}(ii-2) | dsplit{i}(ii) == dsplit{i}(ii-3) | dsplit{i}(ii) == dsplit{i}(ii-4)| dsplit{i}(ii) == dsplit{i}(ii-5)
%                break  
%             elseif ii > 6 && dsplit{i}(ii) == dsplit{i}(ii-1) | dsplit{i}(ii) == dsplit{i}(ii-2) | dsplit{i}(ii) == dsplit{i}(ii-3) | dsplit{i}(ii) == dsplit{i}(ii-4)| dsplit{i}(ii) == dsplit{i}(ii-5) | dsplit{i}(ii) == dsplit{i}(ii-6)
%                 break  
%             elseif ii > 7 && dsplit{i}(ii) == dsplit{i}(ii-1) | dsplit{i}(ii) == dsplit{i}(ii-2) | dsplit{i}(ii) == dsplit{i}(ii-3) | dsplit{i}(ii) == dsplit{i}(ii-4)| dsplit{i}(ii) == dsplit{i}(ii-5) | dsplit{i}(ii) == dsplit{i}(ii-6) | dsplit{i}(ii) == dsplit{i}(ii-7)
%                 break  
              elseif ii == length(dsplit{i})
                 victory = 1;
            end
        end
     end
end
%%% CONCAT all matrices and output
condition_matrix = [];
for i = 1:com_factor
    condition_matrix = [condition_matrix; dsplit{i}];
end
return


function condition_matrix=counter(g, points)
% Function version of old counter.m script.
% No idea who wrote it. Don't fully understand it.
% I assume there's a limit in number of trials it can
% counterbalance, but I have it work with n > 10
% though slow as hell. -DDW Mar/11
% Counter.m
% Creates a counterbalanced sequence of 'g' different trial types
% and 'points' trials.
% usage:  counter(g,points) or [s,k,k2,k3]=counter(g,points)
%         e.g. a = counter(3,60)
%         This will generate a list of 60 items
%              (3 trial types, 20 items per trial type)
%
% s is the counterbalanced sequence of 'points' trials
% k is a vector of the number of each trial type
% k2 is a 2-dimensional matrix of the number of 2-trial combinations
% k3 is a 3-dimensional matrix of the number of 3-trial combinations
% different from counter4 in that g=2 sequences have added randomness
R=5;            % R: added randomness for g=2 sequences (R=1 Random, R>1 Less Random)
                % R can only take on positive integer values
f=1;            % a flag to do the while loop at least once
k=zeros(g);     % because they need to exist before we clear them.
k2=zeros(g,g);         
k3=zeros(g,g,g);
k4=zeros(g,g,g,g);
%%%Send warning
if g > 8
fprintf(['\nWARNING: You are requesting counterbalancing on %d conditions.\n' ...
         'Counterbalancing is offensively slow with 9 or more conditions.\n'], g); 
end
%%%Loop until best possible counter-balanced sequence. Necessary because
%%%of the stochatic nature of the routine.
%%%-Choose the first 'while' statement for quickest execution
%%%-Choose the second 'while' statement for a better counterbalanced
%%% sequence but longer execution.
%%%-Choose the third 'while' statement for the best, but be warned, it
%%% might take a while, a long while.

while(((max(max(k2))-min(min(k2)))>1)|((max(k)-min(k))>1)|f==1);
%while(((max(max(k2))-min(min(k2)))>1)|((max(k)-min(k))>1)|((max(max(max(k3)))-min(min(min(k3))))>1)|f==1);
%while(((max(max(k2))-min(min(k2)))>1)|((max(k)-min(k))>1)|((max(max(max(k3)))-min(min(min(k3))))>1)|(max(max(max(max(k4))))-min(min(min(min(k4)))))>1|f==1);
    %%%Clear and Zero some variables
    clear s k k2 k3 k4 j;
    f=0; k=zeros(1,g); k2=zeros(g,g); k3=zeros(g,g,g); k4=zeros(g,g,g,g);
    %%%Start the crazy
    s(1)=rem(round(rand*10000),g)+1;     %Pick 1st 2 points at random
    s(2)=rem(round(rand*10000),g)+1;
    i=2;
    k(s(1))=k(s(1))+1;k(s(2))=k(s(2))+1;
    k2(s(i-1),s(i))=k2(s(i-1),s(i))+1;   %increment combo counter
    [min2,ind]=min(k2(s(i),:));          % Find the minimum 2-trial combo beginning with s(i)
    same2=find(k2(s(i),:)==min2);        % Check is that minimum is unique
    %%%Check if unique, if so pick it
    if length(same2)==1                 
        s(i+1)=ind;                      % ind is the indice of the minimum.
    else
        s(i+1)=same2(rem(round(rand*10000),length(same2))+1);  % If not unique, pick it at random.
    end
    k(s(i+1))=k(s(i+1))+1;
    k2(s(i),s(i+1))=k2(s(i),s(i+1))+1;   % increment combo counter
    k3(s(i),s(i+1),s(i-1))=k3(s(i),s(i+1),s(i-1))+1;
    %%%The main Loop...
    for i=3:points-1
        [min2,ind]=min(k2(s(i),:));               % same set of statements as above 2-trial type sequences.
        same2=find(k2(s(i),:)==min2);
        if rem(i,R) == 0 & g==2                   % Added 'randomness' for
            s(i+1)=rem(round(rand*10000), g)+1;   % Sequences converge to a periodic signal w/o it.
        elseif length(same2)==1
            s(i+1)=ind;
        else                                          % If not unique, look at 3-trial combos
            [min3,ind]=min(k3(s(i),same2,s(i-1)));    % Find minimum 3-trial combo s(i-1),s(i),same2
            same3=find(k3(s(i),same2,s(i-1))==min3);  % Check if that minimum is unique
            if length(same3)==1                       % If unique, pick it
                s(i+1)=same2(ind);                    % confusing dereferencing...but correct
            else                                      % If not, look at 4-trial combos
                [min4,ind]=min(k4(s(i),same3,s(i-1),s(i-2)));    
                %Min[s(i-2),s(i-1),s(i),same3]
                same4=find(k4(s(i),same3,s(i-1),s(i-2))==min4);  % Unique or not
                if length(same4)==1
                    s(i+1)=same2(same3(ind));                    % more dereferencing...
                else                                             % Pick at Random if not unique
                    s(i+1)=same2(same3(same4(rem(round(rand*10000),length(same4))+1)));
                end
            end
        end
        k(s(i+1))=k(s(i+1))+1;
        k2(s(i),s(i+1))=k2(s(i),s(i+1))+1;
        k3(s(i),s(i+1),s(i-1))=k3(s(i),s(i+1),s(i-1))+1;
        k4(s(i),s(i+1),s(i-1),s(i-2))=k4(s(i),s(i+1),s(i-1),s(i-2))+1;
    end;
    k2;
end;
%Combo_Matrix=k2;     %never used
condition_matrix=s;   %assign output to condition matrix.
return

