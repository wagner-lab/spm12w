function roidata = spm12w_roitool(varargin)
% spm12w_roitool('roi_file','sids','coords')
%
% Inputs
% ------
% roi_file: File specifying the parameters for roi analysis and specification
%           (e.g., 'roi_tutorial.m'). If the path is left unspecified,
%           spm12w_roitool will look in the scripts directory.
%
% sids:     A cell array of Subject IDs for roi analysis. If left unspecified
%           a dialog box will appear asking the user to select subjects. Use 
%           the keyword 'allsids' and all subjects in the specified 
%           glm directory will be used for roi analysis. <optional>
%
% spec:     Cell array of coordinates and roi sizes ([x,y,z,size in mm]) for 
%           manually specified ROI. These will replace any roi specifications 
%           in the roi parameters file. <optional> 
%
% Returns
% -------
% roidata:  A structure containing the fields xxxx. 
%
% spm12w_roitool will generate spherical rois or use pre-existing img masks
% to extract parameter estimates from previously generated contrast files during
% 1st level glm analysis. These parameter estimates may then be submitted to
% a number of simple statistcal tests. These are: 
%       - descriptives: means, standard deviations, min/max values.
%       - ttest1: one-sample t-test across participants
%       - ttest2: idependent sample t-test across participants
%       - correl1: correlation between parameter estimates and subject variables
%       - correl2: same as correl1 but split by group.
%
% Variables for each sid are required for ttest2 and correl1 and correl2
% and are to be specified in a variable file references in the roi
% parameters file (i.e., roi.var_file).
%
% The first argument is the name of a roi parameters file (e.g., roi_tutorial.m).
% The second argument (optional) is a cell array of sids (can be left blank to 
% manually select them). The third argument allows the user to manually 
% specify an roi thereby overriding the rois specified in the roi paramters
% file. Parameter estimates for each specified contrast will be saved to a
% comma seperated csv file formatted for importing into any offline
% statistical software (e.g., R, spss, etc.). A seperate text file
% containing results of basic statistics will also be saved to the roi dir.
%
% Note: If you choose to use the sphere rois genereated by spm12w_roitool
% be aware that if these do not fit the voxel grid, they will be interpolated
% according to interpolation number in the resample variable (default value
% is 0 which is equal to nearest neighbor).
%
% Examples:
%
%       >>spm12w_roitool
%       >>spm12w_roitool('roi_file', './scripts/username/roi_tutorial.m', ...
%                        'sids', {'allsids'}, ...
%                        'spec',{[30,30,21,8]; [22,22,19,6]})
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: January, 2010 | Updated: April, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('roi_file','','sids','','spec','');
args = spm12w_args('nargs',0, 'defaults', args_defaults, 'arguments', varargin);

% Load roi parameters
roi = spm12w_getp('type','roi', 'para_file',args.roi_file);

% Check for cell in case user provided allsids as string.
if ~iscell(args.sids) && ~isempty(args.sids)
    args.sids = cellstr(args.sids);
end

% If sids argument was not provided, open dialog window to get sids.
% If sids argument contained the keyword 'allsids', then get all sids.
% Since we should only do rfx on computed glms, let's look in rfx.glmdir.
if isempty(args.sids)
    args.sids = spm12w_getsid(roi.glmdir);
elseif numel(args.sids) == 1 && strcmp(args.sids,'allsids')
    sids = dir(fullfile(roi.glmdir,'s*'));
    args.sids = {sids.name};
end

% Setup directories for roi analysis. 
spm12w_dirsetup('dirtype','roi','params',roi);

% Parse roi specs, override specs in roi file if custom specs provided
if isempty(args.spec)
    % Check for custom spec csv file and append to roi structure if it exists.
    if exist(fullfile(roi.roispecdir,roi.spec_file),'file')
        spm12w_logger('msg',sprintf('Loading roi specs from csv file: %s',...
            roi.spec_file),'level',roi.loglevel)
        % read in csv as mixed data cell (can't use xlsread on linux so ...)
        csvspec = spm12w_readcsv('csvfile',fullfile(roi.roispecdir,roi.spec_file));
        
    end
    % Parse the csv file
    roistruct = struct();
    for csv_row = 2:size(csvspec,1)
        fieldname = csvspec{csv_row,1};
        if isnumeric(csvspec{csv_row,2}) && ~isempty(csvspec{csv_row,2})
            spec = [csvspec{csv_row,2:end}];
        else
            % Must be nifti mask string, find column containing string
            stridx = find(cellfun(@ischar,csvspec(csv_row,2:end)));
            spec = [csvspec{csv_row,stridx+1}];
        end        
        roistruct.(fieldname) = spec;
    end
    % Append to existing roi specs. If there are none, make it. 
    if isfield(roi,'roi')
        spm12w_logger('msg',['[DEBUG] Merging roi specs from csv ',...
            'and roi parameters files.'],'level',roi.loglevel)
        % add fields of roistruct to roi.roi. Matlab has no direct meanss of
        % merging structures, so this requires a little effort.
        roifieldnames = [fieldnames(roi.roi); fieldnames(roistruct)];
        roi.roi = cell2struct([struct2cell(roi.roi); struct2cell(roistruct)],...
                              roifieldnames, 1);      
    else
        roi.roi = roistruct;
    end  
else
    spm12w_logger('msg',sprintf(['Custom roi specs provided, overriding ',...
              'roi specs in file: %s'],spm_str_manip(roi.para_file,'t')),...
              'level',roi.loglevel)
    % Replace roi spec with the provided coordinates.
    if ~iscell(args.spec)
        args.spec = {args.spec}; % in case user provided single coord
    end
    roistruct = struct();
    for spec = args.spec
        fieldname = sprintf('region_%d_%d_%d',spec{1}(1:3));
        roistruct.(fieldname) = spec{1};
    end
    roi.roi = roistruct;
end

% Check that all appropriate glm dirs and contrasts exist and build up a
% cell array of files per subject per glm condition. Technically this 
% need not be per subject, but just in case different subjects got different
% glms somehow, we should figure out the appropriate con file on a subject
% by subject basis. 
roi.roidata = {}; % Init var for confilenames associated with roi.conds
% Verify the glms and build up the roi.roidata variable.
spm12w_logger('msg',['[DEBUG] Verifying that glm and contrasts exist for ',...
                    'each subject.'],'level',roi.loglevel)
for sid = args.sids
    if ~exist(fullfile(roi.glmdir,sid{1},'SPM.mat'),'file')           
        spm12w_logger('msg',sprintf(['[EXCEPTION] glm ''%s'' for subject %s ', ...
                      'does not exist or is not estimated.'],roi.glm_name, ...
                      sid{1}),'level',roi.loglevel);
        error(['Glm ''%s'' for subject %s does not exist or is not estimated.',...
              'Aborting...'],roi.glm_name,sid{1})      
    else
        for cond = roi.conds
            % Load the SPM file for the GLM (do this subjectwise)
            SPM_ = load(fullfile(roi.glmdir,sid{1},'SPM.mat'));
            % Find the index of the rfx contrast in the SPM.xCon
            conidx = find(strcmp({SPM_.SPM.xCon.name},cond{1}));
            if ~isempty(conidx)
                % Use the index to get the filename of the con file for that sid
                roi.roidata{end+1,1} = sid{1};
                roi.roidata{end,2} = cond{1};
                roi.roidata{end,3} = fullfile(roi.glmdir,sid{1},...
                                              SPM_.SPM.xCon(conidx).Vcon.fname);
            else
                spm12w_logger('msg',sprintf(['[EXCEPTION] condition ''%s'' ',...
                      'for subject %s is not a part of the glm ''%s''.'],...
                      cond{1},sid{1},roi.glm_name),'level',roi.loglevel);
                error(['Condition ''%s'' for subject %s is not part of ',...
                       'glm ''%s''. Aborting...'],cond{1},sid{1},roi.glm_name);      
            end
        end
    end
end

% Resort the file order by contrast and init an empty matrix to hold betas
roi.roidata = sortrows(roi.roidata,[2,1]);

% Iterate through ROIs, pulling data from each contrast and con file.
for roifields = fields(roi.roi)'
    roi_spec = roi.roi.(roifields{1});
    % Check if ROI is a string (i.e., an img mask) or is coordinates
    if ischar(roi_spec)
        spm12w_logger('msg',sprintf(['Performing ROI analysis on region ',...
                      '''%s'' using the mask at:%s'], roifields{1}, roi_spec),...
                      'level',roi.loglevel)
        if ~exist(fullfile(roi_spec),'file')
            if exist(fullfile(roi.roimaskdir, roi_spec),'file')
                roi_spec = fullfile(roi.roimaskdir, roi_spec);
            else
                spm12w_logger('msg',sprintf(['[EXCEPTION] Roi mask ''%s'' ',...
                      'not found. Did you provide the full path?'],...
                      roi_spec),'level',roi.loglevel);
                error('Roi mask %s not found', roi_spec);   
            end        
        end          
        tmpdata = spm12w_readnii('niifiles',roi.roidata(:,3),...
                           'mask',roi_spec,'vox_avg',1,'resample',roi.resample);  
        roi.roidata = [roi.roidata, {tmpdata.maskdata}'];
    elseif isnumeric(roi_spec)
        % Check if roi sphere size is set, if not use default value. 
        if length(roi_spec) ~= 4
            roi_spec(4) = roi.roi_size;
        end
        spm12w_logger('msg',sprintf(['Performing ROI analysis on region ',...
                      '''%s'' using a %dmm sphere at %d,%d,%d'], roifields{1}, ...
                      roi_spec(4), roi_spec(1:3)),'level',roi.loglevel)
        tmpdata = spm12w_readnii('niifiles',roi.roidata(:,3),'sphere',...
                           roi_spec,'vox_avg',1,'resample',roi.resample); 
        roi.roidata = [roi.roidata, {tmpdata.spheredata}'];
    end
end

% Write roidata to tab delimited text file. 
roifile = fullfile(roi.roidir,[roi.roi_name,'.csv']);                               
fid = fopen(roifile,'w'); % Open the file, overwriting prior contents  

% Print the first column headers
fprintf(fid,'%s,%s,%s','Subject','Condition','File');  

% Loop through all ROI regions and print their names as column headers
for region = fields(roi.roi)'
    fprintf(fid,',%s',region{1}); 
end

% Loop through each row of roidata and print subject, condition and file
for roidata_i = 1:size(roi.roidata,1)
    fprintf(fid,'\n%s,%s,%s',roi.roidata{roidata_i,1:3});
    % Now loop through the columns of extracted data
    for roidata_col = 4:size(roi.roidata,2)
        fprintf(fid,',%6.3f',roi.roidata{roidata_i,roidata_col});
    end
end

% Close the file
fclose(fid);
spm12w_logger('msg',sprintf('Roi data saved to roi file ''%s'' at: %s',...
                     [roi.roi_name,'.csv'], roifile),'level',roi.loglevel)
                 
% Generate basic ROI statistics
spm12w_logger('msg',roi.niceline,'level',roi.loglevel)
spm12w_logger('msg',sprintf('Computing basic statistics on data in: %s',...
              [roi.roi_name,'.txt']),'level',roi.loglevel);

% Load and parse the vars file if it exists
if exist(fullfile(roi.roispecdir,roi.var_file),'file')
    spm12w_logger('msg',sprintf('Loading roi variables from csv file: %s',...
        roi.spec_file),'level',roi.loglevel)
    % read in csv as mixed data cell
    csvvar = spm12w_readcsv('csvfile',fullfile(roi.roispecdir,roi.var_file));
    % Remove the first row header to make life easier
    csvhdr = csvvar(1,:);
    csvvar = csvvar(2:end,:);   
    % Get sids order from roidata in case order is different than in var file
    tmpconds = unique(roi.roidata(:,2));
    tmpsids = roi.roidata(ismember(roi.roidata(:,2),tmpconds{1}),1);
    [~,sididx] = ismember(tmpsids,csvvar(:,1));
    % Get group names
    group_names = unique(csvvar(:,2))';
    % Figure out possible numnber of groups (for later)
    comb_grps = nchoosek(1:numel(group_names),2);
    % Search csvvar based on tmpsids and assign groupings
    tmpgrps = csvvar(sididx,2); % groups are now in roidata order
    % Create logical group idx for each group to use for filtering data later
    grp_idx = logical([]);
    for grp = group_names
        grp_idx = [grp_idx, ismember(tmpgrps, grp{1})];
    end
    % Now create correlation variables using sididx to get order right
    var_names = csvhdr(1,3:end);
    var_vars = cell2mat(csvvar(sididx,3:end));
end          

% Set output file for saving statistics.
statsfile = fullfile(roi.roidir,[roi.roi_name,'_stats.txt']);   
diary(statsfile)

% Because we allow on the fly comparisons of conditions, we need to generate a
% list of conditions and interpret any formulas like we do for glm contrasts.
roiconds = {};
for stats = fields(roi.stats)'
    if any(strcmp(roi.stats.(stats{1}),'all_conditions'))
        roiconds = [roiconds, roi.conds];
    else
        roiconds = [roiconds, roi.stats.(stats{1})];    
    end
end
roiconds = sort(unique(roiconds));

% Iterate through these conditions then regions then stat fields
for conds = roiconds
    fprintf('Condition: %s\n', conds{1});
    for region = fields(roi.roi)'
        fprintf('Region: %s\n', region{1});
            % Find which column contains the region data (data starts at col4)
            regioncol = find(ismember(fields(roi.roi), region{1})) + 3;
            % Put the data for region/condition into the data variable.
            % If not a formula, then its simple. 
            if any(ismember(roi.conds,conds{1}))
                % Find the rows containing data for the current condition
                datarows = ismember(roi.roidata(:,2),conds{1});
                % Assign data for this condition and region to a matrix.
                data = cell2mat(roi.roidata(datarows,regioncol));
            else
                % Recycle trick from spm12w_glm_contrast. 
                tmpstr = strsplit(conds{1},'vs.'); %Split on 'vs.' 
                poswt = strtrim(tmpstr{1}); % remove leading/trailing whitespace
                if length(tmpstr) > 1
                    negwt = strtrim(tmpstr{2}); %remove leading/trailing whitespace
                else
                    negwt = {}; % Set to zero if only positives
                end
                % Now split each on wtspace to get members
                poswt = strsplit(poswt);
                % Create pos and neg averages
                posdata = [];           
                for pos_i = 1:numel(poswt)
                    % Find the rows containing data for the current condition
                    datarows = ismember(roi.roidata(:,2),poswt{pos_i});
                    % Assign data for this condition and region to a matrix.
                    posdata(:,pos_i) = cell2mat(roi.roidata(datarows,regioncol));
                end
                % Create neg average if neg elements exist        
                if isempty(negwt)
                    data = mean(posdata,2);
                else
                    negwt = strsplit(negwt);
                    negdata = [];
                    for neg_i = 1:numel(negwt)
                        % Find the rows containing data for the current condition
                        datarows = ismember(roi.roidata(:,2),negwt{neg_i});
                        % Assign data for this condition and region to a matrix.
                        negdata(:,neg_i) = cell2mat(roi.roidata(datarows,regioncol));
                    end
                data = mean(posdata,2) - mean(negdata,2);      
                end                  
            end
        for stats = fields(roi.stats)'
            % Make sure stat is a cell array even if single condition
            stat = cellstr(roi.stats.(stats{1}));
            % Go through stats and check if all_conditions keyword was used
            switch(stats{1})
                case 'descriptives'
                    if any(strcmp(stat, 'all_conditions')) || any(ismember(stat,conds{1}))
                        fprintf('\t=Descriptives:\n');
                        fprintf('\tMean:%4.3f',mean(data));
                        fprintf(' S.D.:%4.3f',std(data));
                        fprintf(' MIN:%4.3f',min(data));
                        fprintf(' MAX:%4.3f\n',max(data));
                    end          
                case 'ttest1'
                    if any(strcmp(stat, 'all_conditions')) || any(ismember(stat,conds{1}))
                        tstats = spm12w_stats('stat','ttest1','y',data);
                        fprintf('\t=one-sample t-test\n');
                        fprintf('\tt-test: t(%d)= %4.2f, p=%4.3f %s\n',...
                                tstats.df,tstats.t,tstats.p,tstats.p_star);
                    end
                case 'ttest2'
                    if any(strcmp(stat, 'all_conditions')) || any(ismember(stat,conds{1}))                       
                        for i_cmb = 1:size(comb_grps,1)     
                            grp1label = group_names{comb_grps(i_cmb,1)};
                            grp2label = group_names{comb_grps(i_cmb,2)};
                            data_grp1 = data(grp_idx(:,comb_grps(i_cmb,1)));
                            data_grp2 = data(grp_idx(:,comb_grps(i_cmb,2)));                            
                            tstats = spm12w_stats('stat','ttest2',...
                                                  'y',data_grp1,'x',data_grp2);
                            fprintf('\t=Independent t-test (%s vs %s)\n',...
                                    grp1label, grp2label);
                            fprintf('\tt(%d)= %4.3f, p=%4.3f %s\n',...
                                    tstats.df,tstats.t,tstats.p,tstats.p_star);
                            fprintf('\tGroup 1 %s: Mean:%4.3f',grp1label, mean(data_grp1));
                            fprintf(' S.D.:%4.3f\n',std(data_grp1));
                            fprintf('\tGroup 2 %s: Mean:%4.3f',grp2label, mean(data_grp2));
                            fprintf(' S.D.:%4.3f\n',std(data_grp2));
   
                        end
                    end
                case 'correl1'
                    if any(strcmp(stat, 'all_conditions')) || any(ismember(stat,conds{1}))
                        for var_i = 1:numel(var_names)
                            rcorr = spm12w_stats('stat','correl','y',...
                                                 data,'x',var_vars(:,var_i));
                            fprintf('\t=Correlation (Var:%s)\n',...
                                    var_names{var_i});
                            fprintf('\tr=%4.2f, p=%4.3f %s\n',...
                                    rcorr.r, rcorr.p, rcorr.p_star); 
                        end
                    end
                case 'correl2'
                    if any(strcmp(stat, 'all_conditions')) || any(ismember(stat,conds{1}))
                        for i_grp = 1:numel(group_names)   
                            data_grp = data(grp_idx(:,i_grp));   
                            if length(data_grp) > 1 % skip correl for n=1
                                for var_i = 1:numel(var_names)
                                    var_grp = var_vars(grp_idx(:,i_grp),var_i);
                                    rcorr = spm12w_stats('stat','correl','y',...
                                                         data_grp,'x',var_grp);
                                    fprintf('\t=Correlation (Var:%s Group:%s)\n',...
                                            var_names{var_i},group_names{i_grp});
                                    fprintf('\tr=%4.2f, p=%4.3f %s\n',...
                                            rcorr.r, rcorr.p, rcorr.p_star); 
                                end            
                            end
                        end
                    end
                otherwise
                    spm12w_logger('msg',sprintf(['[EXCEPTION] Unrecognized ',...
                      'statistics type: %s'],stats{1}),'level',roi.loglevel);
                    error('Unrecognized statistics type: %s.',stats{1});                
            end                   
        end
    end
    fprintf('%s\n',roi.niceline);
end

% Stats complete... Turn off the diary.
diary off

% Save roi mat file
matfile = fullfile(roi.roidir,[roi.roi_name,'.mat']);
save(matfile,'roi');

% Print final words
msglist{1} = roi.niceline;
msglist{2} = sprintf('Roi data   : %s', roifile);
msglist{3} = sprintf('Statistics : %s', statsfile);
msglist{4} = sprintf('Parameters : %s', matfile);

for msg = msglist
    spm12w_logger('msg',msg{1},'level',roi.loglevel)
end    