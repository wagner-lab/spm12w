function csvcell = spm12w_readcsv(varargin)
% spm12w_readcsv('csvfile')
%
% Input
% -----
% csvfile : Full path to csv file.
%
% Returns
% -------
% csvcell : A cell array containing the data in the csvfile. Numeric data will 
%           be converted to numeric cells. 
%
% Simple tool to read in csv files and output a cell array with mixed data types
% (numeric, string). This replicates the function of matlab's xlsread which
% will output a cell array of mixed data types on machines with MS office, but
% when run in basic mode on a linux machine, this "feature" of xlsread no 
% longer functions. So we roll our own!
%
%
% Examples:
%
% Read in a csv spec file
%   >> csvcell = spm12w_readcsv('csvfile','tutorial_spec.csv');
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: May, 2015 | Updated: May, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('csvfile','');
args = spm12w_args('nargs',2, 'defaults', args_defaults, 'arguments', varargin);

% Open to csv file and read every line as strings
fid = fopen(args.csvfile,'r');
csvtmp = textscan(fid,'%s','delimiter','\n');
fclose(fid);

% Parse each row of strings and convert to numeric if necesasry
csvcell = {};
csvrow = 1;
for csv = csvtmp{1}'
    csvcelltmp = regexp(csv{1},',','split');
    for i = 1:length(csvcelltmp)
        [csvnum, numstatus] = str2num(csvcelltmp{i});
        if ~isempty(csvcelltmp{i})
            if numstatus == 1
                csvcell{csvrow,i} = csvnum;
            else
                csvcell{csvrow,i} = csvcelltmp{i};              
            end
        end
    end
    csvrow = csvrow + 1;
end
