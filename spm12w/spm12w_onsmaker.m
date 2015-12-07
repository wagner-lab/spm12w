function parameters = spm12w_onsmaker(varargin)
% spm12w_onsmaker(csvfile)
%
% Input
% -----
% csvfile: csv file with 1 row of headers followed by onsets to export.
%
% Takes a csv formatted file with subjects, conditions, parametrics, 
% TRs and RTs and outputs condition onsets, parametric files, and 
% subject_condition_dur.txt files. Will also converts any NR to mean of 
% condition. To use you must format file accorting to: 
% tutorial_data_onsets.csv.
%
% Onsets will be output a directory named after your csv file. This 
% directory can be found in the same directory as the csv file you 
% provided. You may move and rename this directory afterwards.
%
% Examples:
%
% Without arguments, spm12w will prompt the user for a csv file.
%
%   >> spm12w_onsmaker
%
%   >> spm12w_onsmaker('csvfile','./auxil/tutorial_data_onsets.csv')
%
% TODO: This is a bit of a rush job. Optimize code, improve logging add
% option to output sid specific condition and parametric onsets based on
% head nam tokens. (sid_parname, sid_condname) etc.
% 
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: June, 2010 | Updated: December, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
arg_defaults = struct('csvfile','');
args = spm12w_args('nargs',0, 'defaults',arg_defaults, 'arguments',varargin);

% Select csv file if not provided then move to csv file directory and load
if isempty(args.csvfile)
    msg = 'Please select your csv file containing onsets.';
    args.csvfile = spm_select(1,['','','_.*\.csv$'],msg,[]);
end

% in case user supplies a file in current dir.
if isdir(fullfile(spm_str_manip(args.csvfile,'h')))
    cd(fullfile(spm_str_manip(args.csvfile,'h')));
end

% Open csv file and parse to a num and txt format.
fid = fopen(args.csvfile, 'r');
    header = fgetl(fid);
    % parse header
    i = 1;
    while ~isempty(header)
        [t, header] = strtok(header,',');
        txt(i) = {t};
        i = i + 1;
    end
    % parse for txt
    col2 = textscan(fid, '%*s%s%*[^\n]','Delimiter',',');
    % pase for num
    frewind(fid);
    fgetl(fid); %skips header
    numtmp = textscan(fid, ['%n%*s',repmat('%n',1,length(txt)-2)],'Delimiter',',','TreatAsEmpty',{'NA','na','NR','nr'});
fclose(fid);
clear ans fid header i onsmatrix t 
% Reconstruct txt and num to match the old xlsread output
% Put col2 into column 2 of txt var
txt(2:length(col2{1})+1,2) = col2{1}(1:end);
% Fix num so that it includes NANs in 2nd column
num      = numtmp{:,1};
num(:,2) = NaN(length(num),1);
for i = 2:length(numtmp)
    num(:,i+1) = numtmp{:,i};
end

% Parse txt data and assign names and condition sizes to structure.
% Put subject names in an array
fprintf('Submitting onset data to a rather invasive search...');

parstart = '';
sidstart = '';
% Get parametric starting location (if any)
for parloc = 1:size(txt,2)
    % tr at 1, cond at 2, par or sid specific could start at 3. 
    if regexp(txt{1,parloc},'par_')
        parstart = parloc;
        break    
    end
end

% Assuming onsfile is formatted correctly, sid should start after
% parametrics and be of format s##. Iterate through txt looking for sids
for sidloc = parstart:size(txt,2)
    % tr at 1, cond at 2, sid specific could start at 3. search for sid
    % format. Search using reg expression.
    if regexp(txt{1,sidloc},'s[0123456789][0123456789]')
        sidstart = sidloc;
        break    
    end
end

% Assign parnames (if exist) and sids (if exist) to onsdata structure
if ~isempty(parstart) && ~isempty(sidstart)
    % clean up parnames
    onsdata.parnames = strrep(txt(1,parstart:sidstart-1),'par_','x');
    onsdata.sids = txt(1,sidstart:end);
elseif ~isempty(sidstart)
    onsdata.sids = txt(1,sidstart:end);
end

% Put conditions and condition length into an array (should be col2)
tmp = txt(2:size(txt,1),2);
condname(1) = tmp(1);
condlength(1) = 1;
current=1;
for i = 1:length(tmp)
    if i > 1
        if strcmp(condname{current},tmp{i})
            condlength(current)= condlength(current) + 1;
        else
            current=current+1;
            condname(current) = tmp(i);
            condlength(current) = 1;
        end
    end
    fprintf('.'); 
end
% dump to structure
onsdata.cond=condname;
onsdata.condsize=condlength;

% Assign TRs and parametrics to appropriate conditions
iTR = 1;
for icond = 1:length(onsdata.cond)
    onsdata.tr{icond}=num(iTR:iTR+onsdata.condsize(icond)-1,1);
    parcount = 1;
    for ipar = parstart:sidstart-1
        onsdata.par{icond, parcount} = num(iTR:iTR+onsdata.condsize(icond)-1,ipar);
        parcount = parcount + 1;
    end
    iTR=iTR+onsdata.condsize(icond);
    fprintf('.');
end

% Correct RTs for NANs and assign to structure
fprintf('\nReplacing missing RTs with within condition mean');
% Make a data structure
iTR=1;
for i = 1:length(onsdata.cond)
    onsdata.rtdata{1,i}=num(iTR:iTR+onsdata.condsize(i)-1,sidstart:end);    
    %calc mean without nans
    tmpmean=nanmean(onsdata.rtdata{i},1);
    %find out where is nan
    notnan=isnan(onsdata.rtdata{i});
    %replace with mean per subject
    tmpdata=onsdata.rtdata{i};
    for ii = 1:length(onsdata.sids)
        tmpdata(notnan(:,ii),ii)=tmpmean(ii);
        fprintf('.')
    end
    onsdata.rtdata{i}= tmpdata; 
    iTR=iTR+onsdata.condsize(i);
    fprintf('.')
end
fprintf('\nNo responses destroyed and replaced with within condition mean\n');

% write onsets
mkdir(fullfile(spm_str_manip(args.csvfile,'tr')))
cd(fullfile(spm_str_manip(args.csvfile,'tr')))

fprintf('Writing out files\n');
% condition onsets
for i = 1:length(onsdata.cond)
    fname=[onsdata.cond{i},'.txt'];
    tmpdata=onsdata.tr{i};
    save(fname,'tmpdata','-ASCII');
    fprintf('%s written\n',fname);
end

%Parametric onsets 
for ipar = 1:length(onsdata.parnames)
    for icond = 1:length(onsdata.cond)
        %check first row for NAN, if nan, assume all rows are nan and skip
        if ~isnan(onsdata.par{icond,ipar}(1))
            fname=[onsdata.cond{icond},onsdata.parnames{ipar},'.txt'];
            tmpdata=onsdata.par{icond,ipar};
            save(fname,'tmpdata','-ASCII');   
            fprintf('%s written\n',fname);
        end
    end
end

% Duration onsets.
for i = 1:length(onsdata.sids)
    for ii = 1:length(onsdata.cond)
        fname=[onsdata.sids{i},'_',onsdata.cond{ii},'_dur.txt'];
        tmpdata=onsdata.rtdata{ii};
        tmpdata=tmpdata(:,i);
        save(fname,'tmpdata','-ASCII');
        fprintf('%s written\n',fname);
    end
end

fprintf('===onset files written to: %s\n',pwd);
cd('..')
