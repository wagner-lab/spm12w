function spm12w_preprocess(varargin)
% spm12w_preprocess('sub_id','parameter_file')
%
% Inputs
% ------
% sid:            Subject ID of subject to be preprocessed (e.g., 's01')
%
% parameter_file: File specifying the parameters for preprocessing (e.g.,
%                 'p_tutorial.m'). If the full path is left unspecified, 
%                 spm12w_preprocess will look in the scripts directory.
%                 <optional>
%
% spm12w_preprocess will perform image preprocessing for anatomical and
% functional imaging data. The first argument is a sid. The second 
% argument is the name of a parameter files (e.g., p_studyname.m) and is 
% optional. If the parameter file is unspecified, matlab will prompt the
% user to select the file.
%
% Preprocessed data will be output to the prep/preproc_name/sid
% directory. 
%
% A log file (preprocess.log) will be written to the subject's preproc 
% directory and contains information and debug messages from preprocessing 
%
% A mat file (preprocess.mat) will be written to the subject's preproc
% directory containing the structure with all preprocessing parameters
%
% A pdf file (preprocess.pdf) will also be written to the subject's preproc
% directory and contains the matlab figure outputs of individual 
% preprocessing steps.
%
% Examples:
%
%       >>spm12w_preprocess('sid','s01')
%       >>spm12w_preprocess('sid','s01', ...
%                           'para_file','./scripts/username/p_tutorial.m')
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: March, 2006 | Updated: June, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('sid','', 'para_file','');
args = spm12w_args('nargs',2, 'defaults', args_defaults, 'arguments', varargin);

% Load parameters
p = spm12w_getp('type','p', 'sid',args.sid, 'para_file',args.para_file);

% Clean up prepro directory
spm12w_dirsetup('dirtype', 'prep_clean', 'params',p)

% Setup logfile
spm12w_logger('msg','setup_prep', 'level',p.loglevel, 'params',p)

% Get parameters from epi files in rawdir
tmp_epifiles = dir([fullfile(p.rawdir,p.fmri),'_*.nii.gz']);
epifiles = {};
for epifile = {tmp_epifiles.name}
    epifiles = [epifiles, fullfile(p.rawdir,epifile{1})];
end
p = spm12w_getscanner('epifiles',epifiles, 'p',p);

% Setup repro directory with fresh files.
spm12w_dirsetup('dirtype', 'prep_setup', 'params',p)

% cd to data directory 
% Since internal spm_print calls will print figure to cwd, might as well
% be in p.datadir.
cd(p.datadir)

% Open (and hide) an spm graphics figure or else spm_realign won't find it.
F = spm_figure('CreateWin','Graphics', 'spm12w preprocessing', 'off'); 
 
% Begin Preprocessing
spm12w_logger('msg',p.niceline, 'level',p.loglevel)
spm12w_logger('msg',sprintf('Beginning preprocessing on subject: %s', ...
                             p.sid), 'level',p.loglevel)

% Step: Trim volumes
%       Trim volumes from beginning and end of epi files (e.g., to remove 
%       overlapping sections of movie data). p.trimvols can be a scalar(5),
%       vector [5,5] or matrix where rows = nses (e.g., [5,5;10,10]). In
%       the case of scalar or vectors, trim volumes will convert to 
%       nses x 2 matrix. 
if sum(any(p.trimvols)) % check for nonzero
    spm12w_logger('msg',p.niceline, 'level',p.loglevel)
    spm12w_logger('msg',sprintf('Trimming volumes for subject: %s', ...
                  p.sid),'level',p.loglevel)
    % Adjust p.trimvols if user gave single number or single row or colum
    rep_adjust = size(zeros(p.nses,2)) + 1 - size(p.trimvols); %repmat factor
    trimvols = repmat(p.trimvols, rep_adjust); 
    if any(size(zeros(p.nses,2)) ~= size(trimvols))
        spm12w_logger('msg',['[EXCEPTION] Trim volume matrix does not ' ...
                      'match expected size '], 'level',p.loglevel)
        diary off
        error(['Trim volume matrix does not match expected size ' ...
              'based on number of runs...'])    
    end   
    % Trim each session by trimvols .
    for ses = 1:p.nses       
        spm12w_logger('msg',sprintf(['[DEBUG] Triming first %d volumes & last ' ...
                              '%d volumes for run %d of %d'], ...
                              trimvols(ses,1), trimvols(ses,2), ses, ...
                              p.nses),'level',p.loglevel)
        epifile = sprintf('%s_r%02d.nii',p.fmri,ses);
        epitrimfile = sprintf('t%s_r%02d.nii',p.fmri,ses);
        % load hdr and data for range
        vdata = spm12w_readnii('niifile',fullfile(p.datadir,epifile), ...
                               'range', ...
                               [trimvols(ses,1)+1:p.nvols(ses)-trimvols(ses,2)]);
        % write new file of trimmed data
        spm12w_writenii('niifile',fullfile(p.datadir,epitrimfile), ...
                        'hdr',vdata.hdr,'data',vdata.data);                        
        spm12w_logger('msg',sprintf(['[DEBUG] Adjusting nvols to reflect new ' ...
                      'volume length of: %d vols'],length(vdata.hdr)), ...
                      'level',p.loglevel)
        p.nvols(ses) = length(vdata.hdr); % adjust number of TRs
    end
    spm12w_logger('msg','Trimming complete...', 'level',p.loglevel)
    p.fmri = ['t',p.fmri];
end

% Step: Shuffle check 
%       Based on shufflecheck.m by Petr Jananta
if p.shuffle
    spm12w_logger('msg',p.niceline, 'level',p.loglevel)
    spm12w_logger('msg',sprintf('Shuffle check on subject: %s', p.sid), ...
                  'level',p.loglevel)
    for ses = 1:p.nses
        spm12w_logger('msg',sprintf(['[DEBUG] Computing shuffle check on run %d ' ... 
                               'of %d'], ses, p.nses), 'level',p.loglevel)
        epifile = sprintf('%s_r%02d.nii',p.fmri,ses);
        epiV = spm_vol(fullfile(p.datadir,epifile));
        p.shuff{ses} = zeros(p.nslice, p.nvols(ses));
        for ivol = 1:p.nvols(ses)
            for islice = 1:p.nslice    
                tmp = spm_slice_vol(epiV(ivol),spm_matrix([0 0 islice]), ...
                                    epiV(ivol).dim(1:2),0);
                p.shuff{ses}(islice,ivol) = sum(tmp(:));
            end
        end
    end
    spm12w_logger('msg','Shuffle check complete...', 'level',p.loglevel)
end

% Step: Despike
%       Originally ran epi data through 3dDespike, however we're going to
%       switch to a matlab implementation of 3dDespike's algorithm for
%       winpc platform compatability.
%       Matlab 3dDespike:
%                   -Fit a curve (L1 i.e., least absolute deviation reggression) 
%                    to each voxel timeseries
%                   -Compute the MAD of the difference between curve
%                    and the timeseries
%                   -Estimate the standard deviation 'sigma' of the residuals"
%                    as sqrt(PI/2)*MAD.
%                   -For each voxel value, define s = (value-curve)/sigma
%                   -Values with s > c1 are replaced with a value that yields"
%                    a modified s' = c1+(c2-c1)*tanh((s-c1)/(c2-c1))
%                   -c1 is the threshold value of s for a 'spike' [default
%                   c1=2.5].
%                   -c2 is the upper range of the allowed deviation from the curve:\n"
%                   s=[c1..infinity) is mapped to s'=[c1..c2)   [default c2=4]

if p.despike
    for ses = 1:p.nses
        spm12w_logger('msg',p.niceline, 'level',p.loglevel)
        spm12w_logger('msg',sprintf('Despike on subject: %s', p.sid),...
                     'level',p.loglevel)
        epifile = sprintf('%s_r%02d.nii',p.fmri,ses);
        dsfile  = sprintf('d%s_r%02d.nii',p.fmri,ses);
        spm12w_despike('niifile',fullfile(p.datadir,epifile), ...
                       'outfile',fullfile(p.datadir,dsfile), 'cut', p.cut,...
                       'c_order',p.c_order, 'polort', p.dpolort', ...
                       'least_sq', p.least_sq, 'dilate', p.dilate, ...
                       'loglevel', p.loglevel);
    end
        p.fmri = ['d',p.fmri];
end

% Step: Slice Time Correction
if p.slicetime
    spm12w_logger('msg',p.niceline, 'level',p.loglevel)
    spm12w_logger('msg',sprintf('Slice time correction on subject: %s', p.sid),...
                 'level',p.loglevel)    
    % Slice time correction now allows matching of sessions with different TRs
    % This complicates things as we have to run slicetime correction
    % everytime we detect a new TR. Number of sessions = unique(p.tr).
    for tr = num2cell(unique(p.tr))
        % custom log message for the rare case when unique(p.tr) > 1
        if length(unique(p.tr)) > 1
            spm12w_logger('msg',sprintf(['Slice time correction for sessions ',...
                      'with tr: %.1f'], tr{1}),'level',p.loglevel)
        end
        % Load epi files for files corresponding to the current unique tr      
        epifiles = cell(1,sum(ismember(p.tr,tr{1})));
        epi_i = 1;
        for ses = 1:p.nses
            if p.tr(ses) == tr{1} 
                epifile = sprintf('%s_r%02d.nii',p.fmri,ses);
                epifiles{epi_i} = fullfile(p.datadir,epifile);
                epi_i = epi_i + 1;
            end
        end  
        % Perform slice time correction
        spm12w_logger('msg',sprintf('TR for current sessions: %.1f', tr{1}), ...
                     'level',p.loglevel)  
        spm12w_logger('msg',sprintf('Reference slice: %d', p.refslice), ...
                     'level',p.loglevel)  
        spm12w_logger('msg',sprintf('Slice time acquisition order: %s', ...
                      mat2str(p.sliceorder)), 'level',p.loglevel)                
        TA = tr{1}-tr{1}/p.nslice;    
        spm_slice_timing(epifiles,p.sliceorder,p.refslice, ...
                        [TA/(p.nslice-1) tr{1}-TA]);                     
    end    
    spm12w_logger('msg','Slice time correction complete...', 'level',p.loglevel)
    p.fmri = ['a',p.fmri];
end

% Step: Realign epi images
if p.realign
    spm12w_logger('msg',p.niceline, 'level',p.loglevel)
    spm12w_logger('msg',sprintf('Realignment on subject: %s', p.sid),...
                 'level',p.loglevel)  
    % Load epi files
    epifiles = cell(1,p.nses);
    for ses = 1:p.nses
        epifile = sprintf('%s_r%02d.nii',p.fmri,ses);
        epifiles{ses} = fullfile(p.datadir,epifile);
    end 
    % Set realign flags and run realignment.
    realign_flags = struct('quality',p.realign_quality, 'rtm',p.realign_rtm, ...
                            'wrap',p.wrap_r, 'interp',p.interp_r); 
    spm_realign(epifiles, realign_flags);
    % delete spm's ps file (as it contains current date and will not append if
    % prepro extends into the next day). 
    delete(sprintf('spm_%s.ps', datestr(now,'yyyymmmdd')))
    print(F, 'preprocess.ps', '-dpsc2','-painters','-append','-noui')
    % NB: No need to write our realigned data if using transformation matrix at 
    % the next reslice (i.e. unwarping or normalization). Be careful if
    % running cleanup prior to Dartel as it could harm the realignment if
    % mat files are missing for Dartel (maybe don't delete mat files?). 
    % Create mean image if no unwarping is selected (otherwise unwarp does it). 
    if ~p.unwarp && ~strcmp(p.normalize,'none')
        flags = struct('which',[0 1]); %Sets to no reslice (0) output mean (1)
        spm_reslice(epifiles,flags);
    elseif ~p.unwarp && strcmp(p.normalize,'none')
        flags = struct('which',[1 1]); 
        spm12w_logger('msg',sprintf(['Writing out realignment images for '...
                      'subject: %s'], p.sid), 'level',p.loglevel)  
        spm_reslice(epifiles,flags); 
    end
    spm12w_logger('msg','Realignment complete...', 'level',p.loglevel)
    % Make shufflecheck figure (we make it here because we need rp)
    if (p.shuffle)
        ys = (1:p.nslice)*0.25+4;  %yticks
        % Layout depends on number of sessions (big for 3, small for >3)
        subplot_ses = 1;
        for ses=1:p.nses
            ra = load(fullfile(p.datadir,sprintf('rp_%s_r%02d.txt',p.fmri,ses)));
            p.ra{ses} = ra;
            if p.nses <= 2
                subplot(2,1,subplot_ses); 
            elseif p.nses == 3
                subplot(3,1,subplot_ses); 
            else 
                subplot(3,2,subplot_ses); 
            end
            imagesc(1:p.nvols(ses),ys,p.shuff{ses});
            set(gca,'YDir','Normal');
            hold on
            colormap jet
            dra   = diff(ra(:,1:3));
            [x,~] = find(abs(dra)>1);
            plot(ra-repmat(ra(1,:),p.nvols(ses),1));
            for j=1:length(x)
                plot([x(j),x(j)],[-5,20],'k--');
            end
            axis([0,p.nvols(ses),-4,max(ys)]);
            title(sprintf('Subject: %s epi run: %d', p.sid, ses));
            % If 6 runs on page, print page, clear and start over for
            % next set of 6 runs.
            if ~mod(ses,6)
                print(F, 'preprocess.ps', '-dpsc2','-painters','-append','-noui')
                spm_figure('clear',F)  
                subplot_ses = 1;  %reset the subplot
            else
                subplot_ses = subplot_ses + 1;
            end 
        end
        % Print last shuffle check figure to ps file
        print(F, 'preprocess.ps', '-dpsc2','-painters','-append','-noui')
    end
    % Set bold token to rbold if realigned images were written
    if ~p.unwarp && strcmp(p.normalize,'none')
        p.fmri = ['r',p.fmri];
    end 
end

% Step: Unwarp
if p.unwarp
    spm12w_logger('msg',p.niceline, 'level',p.loglevel)
    spm12w_logger('msg',sprintf('Unwarping subject: %s', p.sid),...
                 'level',p.loglevel)  
    % Load epi files
    epifiles = cell(1,p.nses);
    for ses = 1:p.nses
        epifile = sprintf('%s_r%02d.nii',p.fmri,ses);
        epifiles{ses} = fullfile(p.datadir,epifile);
    end         
    % Estimate unwarping parameters
    close(F) % Or else we get hit with figures stealing focus.
    for ses=1:p.nses
        spm12w_logger('msg',sprintf(['[DEBUG] Estimating unwarping ' ...
                      'parameters for run %d of %d'],ses, p.nses), ...
                      'level',p.loglevel)  
        epitmp              = spm_vol(epifiles{1});
        uwe_flags.M         = epitmp(1).mat;
        ds                  = spm_uw_estimate(epifiles{ses}, uwe_flags);   
        ads(ses)            = ds;
        [path,name]         = fileparts(epifiles{ses});
        pefile              = fullfile(path,[name '_uw.mat']);
        save(pefile,'ds');
    end
    % Write unwarped images (more efficient if spm_uw_apply returned VO)
    spm12w_logger('msg','Writing out unwarped images','level',p.loglevel)  
    uwr_flags = struct(  'wrap',        p.wrap_r,...   %Default is [0 0 0]
                       'interp',      p.interp_r);     %Default is 4  
    spm_uw_apply(ads,uwr_flags);
    spm12w_logger('msg','Unwarping complete...', 'level',p.loglevel)
    p.fmri = ['u',p.fmri]; 
    % Activate F again so as not to interfere with future figures
    F = spm_figure('CreateWin','Graphics', 'spm12w preprocessing', 'off'); 
    spm_figure('Clear',F);
end

%Step: Mode 1000 Normalization (WIP)
if isfield(p,'gms1k') && p.gms1k == 1 || isfield(p,'gms1k') && p.gms1k == 2
    spm12w_logger('msg',p.niceline, 'level',p.loglevel)
    spm12w_logger('msg',sprintf('Mode 1000 normalization for subject: %s', ...
                  p.sid), 'level',p.loglevel)  
    % Load epi files
    epifiles = cell(1,p.nses);
    for ses = 1:p.nses
        epifile = sprintf('%s_r%02d.nii',p.fmri,ses);
        epifiles{ses} = fullfile(p.datadir,epifile);
    end 
    %%% Calculate Mean or Media
    %%% Notes for myself: spm_global calculates the global mean over a 3D
    %%% image. It's in mex/c code unfortunately so can't view m source.
    %%% Reverse engineering it: For a given image, spm_global gives 180.04
    %%% It first calculates the mean (including zeros) then calcualtes
    %%% it again by taking only elements > mean/8. We can replicate this
    %%% function using mean(z1(z1>mean2(z1)/8)) or mean(z1(z1>mean(z1(:))/8) 
    %%% which gives 180.04. So that's what the global does, the mode would 
    %%% work the same way. (by the way this works for 3d or 4d). But be
    %%% careful when claculating variance. var of z1(:) can't use that trick.

    %epitmp = epifiles{1};
    %img = spm_read_vols(spm_vol(epitmp));
    %   
    %   
    %   z1   = img(:,:,:,1)
    %   z2   = img(:,:,:,150)
    %   
    %   
    %   spm_global(z1) = 180.0403
    %   mean of z1(z1~=0) = 109.7267
    %   mean(z1(z1>109.7267/8)) = 211.8839 %too big
    %   mean2(z1) = 51.6278 %includes zeros
    %   mean(z1(z1>51.6278/8)) = 180.0403 %so this is what spm_global does!
    %mean 2 will mean entire matrix
    spm12w_logger('msg','Mode 1000 normalization complete...',...
                  'level',p.loglevel)
    p.fmri = ['m',p.fmri]; 
end

%Step: Normalise EPI images
if ~strcmp(p.normalize, 'none')
    spm12w_logger('msg',p.niceline, 'level',p.loglevel)
    spm12w_logger('msg',sprintf('Normalizing (type: %s), subject: %s', ...
                  p.normalize, p.sid),'level',p.loglevel)     
    % Load fileanmes for last preprocessing stage
    epifiles = cell(1,p.nses);
    for ses = 1:p.nses
        epifile = sprintf('%s_r%02d.nii',p.fmri,ses);
        epifiles{ses} = fullfile(p.datadir,epifile);
    end 
   
    switch(p.normalize)
        case 'epi'
            % Load last preprocessed stage
            V = spm_vol(epifiles);
            V = cat(1,V{:});
            % Memory map the mean image 
            meanf   = fullfile(p.datadir, sprintf('mean%s_r01.nii', p.fmri));
            Vm      = spm_vol(meanf);
            % Normalization
            matname  = [spm_str_manip(Vm.fname,'sd') '_sn.mat'];
            VG = fullfile(spm('Dir'),'toolbox','OldNorm','EPI.nii');
            spm12w_logger('msg','[DEBUG] Determining normalization parameters', ... 
                          'level',p.loglevel)
            % EPI to EPI normalization is depreciated in spm12, so adjust path to add it
            oldnorm = fullfile(spm('Dir'),'toolbox','OldNorm');
            addpath(oldnorm)
            % get normalization parameters using spm's defaults (I see no reason to
            % divert from the old normalise spm defaults). 
            params   = spm_normalise(VG,Vm,matname,'','', ...
                       spm_get_defaults('old.normalise.estimate')); 
            % delete spm's ps file and print append to our own ps file
            delete(sprintf('spm_%s.ps', datestr(now,'yyyymmmdd')))
            print(F, 'preprocess.ps', '-dpsc2','-painters','-append','-noui')       
            % Set normalization defaults, then normalize mean and epi files 
            defaults.normalise.write.vox    = p.evoxsize;
            defaults.normalise.write.bb     = p.boundbox;
            defaults.normalise.write.interp = p.interp_w;
            defaults.normalise.write.wrap   = p.wrap_w;
            spm12w_logger('msg','[DEBUG] Writing out normalized images', ... 
                          'level',p.loglevel)
            msk      = spm_write_sn(V,params,defaults.normalise.write,'mask');
            spm_write_sn(Vm,params,defaults.normalise.write,msk);	% write nmean
            % Write normalised epi images
            spm_write_sn(V,params,defaults.normalise.write,msk);
            % Adjust path back to what it was
            rmpath(oldnorm)    

        case 'spm12'
            if p.coreg2epi % you would be a fool not to, but check anyway.
                spm12w_logger('msg','Coregistering anatomy to epi data...', ... 
                              'level',p.loglevel)                
                % copy the anatomy
                copyfile('anat.nii', 'canat.nii');
                % Setup coregistration
                meanf = fullfile(p.datadir, sprintf('mean%s_r01.nii', p.fmri));
                target = spm_vol(meanf);
                source = spm_vol('canat.nii');
                coreg_flags = spm_get_defaults('coreg.estimate');            
                coreg_x = spm_coreg(target,source,coreg_flags);              
                % Write coregistration transform to sform field in header
                M   = spm_matrix(coreg_x);         %convert coreg trans to matrix
                MMS = spm_get_space(source.fname); %Get current image space
                spm_get_space(source.fname, M\MMS);%Set transform in header.
                anat = 'canat.nii';
                spm12w_logger('msg',['Coregisteration transformation saved ',...
                              'to header of file:canat.nii'],'level',p.loglevel) 
                % delete spm's ps file and print append to our own ps file
                delete(sprintf('spm_%s.ps', datestr(now,'yyyymmmdd')))
                print(F, 'preprocess.ps', '-dpsc2','-painters','-append','-noui')   
            else
                anat = 'anat.nii';
            end
            % Adjust paths to include spm12/config dir
            addconfig = fullfile(fileparts(which('spm.m')),'config');
            addpath(addconfig)
            % Estimate normalization deformations based on structural anat
            spm12w_logger('msg',sprintf(['Estimating normalization ',...
                          'deformations based on file:%s'], anat), ... 
                          'level',p.loglevel) 
            job.subj.vol = {sprintf('%s,1',anat)};
            job.eoptions.biasreg = p.ebiasreg;
            job.eoptions.biasfwhm = p.ebiasfwhm;
            job.eoptions.tpm = {p.etpm};
            job.eoptions.affreg = 'mni';
            job.eoptions.reg = [0 0.001 0.5 0.05 0.2];
            job.eoptions.fwhm = 0;
            job.eoptions.samp = p.esamp;
            def_fname = spm_run_norm(job);
            % Apply deformations to anat at avoxsize
            spm12w_logger('msg',sprintf(['Applying normalization ',...
                  'deformations to anatomy file:%s'], anat), ... 
                  'level',p.loglevel)  
            job.subj.def = def_fname.def; % the estimated deformation name
            job.subj.resample = {sprintf('%s,1',anat)}; % file to normalize
            job.woptions.prefix = 'w'; %new to spm12_6685
            job.woptions.bb = p.boundbox;
            job.woptions.vox = p.avoxsize;    % voxelsize for anatomy
            job.woptions.interp = p.interp_w; % we love those big splines
            spm_run_norm(job);
            % Apply deformations to epi data at avoxsize
            spm12w_logger('msg',['Applying normalization deformations to ',...
                  'epi files'],'level',p.loglevel)
            job.subj.def = def_fname.def;  
            job.subj.resample = epifiles'; %need to transpose
            job.woptions.bb = p.boundbox;
            job.woptions.vox = p.evoxsize;
            job.woptions.interp = p.interp_w;
            spm_run_norm(job);
            % rm spm12/config from path
            rmpath(addconfig)
            
        case 'dartel'
            % do nothing for now.
            
    end
    spm12w_logger('msg',sprintf('Normalizing (type: %s) complete...',...
                  p.normalize),'level',p.loglevel)
    p.fmri = ['w',p.fmri];       
end

% Step: Smoothing
if p.smoothing
    spm12w_logger('msg',p.niceline, 'level',p.loglevel)
    spm12w_logger('msg',sprintf('Smoothing (%dmm FWHM) subject: %s', ...
                  p.smoothing, p.sid),'level',p.loglevel) 
    for ses = 1:p.nses
        epifile = fullfile(p.datadir,sprintf('%s_r%02d.nii',p.fmri,ses));
        smfile = fullfile(p.datadir,sprintf('s%s_r%02d.nii',p.fmri,ses));
        spm_smooth(epifile, smfile, p.smoothing)
    end 
    spm12w_logger('msg',sprintf('Smoothing (%dmm FWHM) complete...', ...
                  p.smoothing), 'level',p.loglevel)
    p.fmri = ['s',p.fmri];
end

% Step: SNR 
if p.snr
    spm12w_logger('msg',p.niceline, 'level',p.loglevel)
    spm12w_logger('msg',sprintf('Calcualating tSNR on subject: %s', ...
                  p.smoothing, p.sid),'level',p.loglevel) 
    % Load last preprocessing stage
    epifiles = cell(1,p.nses);
    for ses = 1:p.nses
        epifile = sprintf('%s_r%02d.nii',p.fmri,ses);
        epifiles{ses} = fullfile(p.datadir,epifile);
    end 
    spm12w_snr('niifiles',epifiles,'outdir',p.qadir,'sid',p.sid,...
               'psname','preprocess.ps','detrend',0,'loglevel',p.loglevel);
    spm12w_logger('msg','tSNR calculation complete...', 'level',p.loglevel) 
end

% Step: Slice Noise Analysis 
%       Borrowed from slices_analyse by Antonia Hamilton
if p.slices
    spm12w_logger('msg',p.niceline, 'level',p.loglevel)
    spm12w_logger('msg',sprintf(['Calcualating slice noise check on ' ...
                  'subject: %s'], p.sid),'level',p.loglevel) 
    epifiles = cell(1,p.nses);
    for ses = 1:p.nses
        if p.smoothing > 0
            epifile = sprintf('%s_r%02d.nii',p.fmri(2:end),ses);
        else
            epifile = sprintf('%s_r%02d.nii',p.fmri,ses);
        end
        epifiles{ses} = fullfile(p.datadir,epifile);
        if ~isfield(p, 'ra')
            p.ra{ses} = '';
        end
    end 
    % Open a new figure (as it gets closed somewhere)
    % or else we get focus stolen by slicenoise
    F = spm_figure('CreateWin','Graphics', 'spm12w preprocessing', 'off'); 
    spm12w_slicenoise('niifiles',epifiles,'radata',p.ra,'mask',p.mask, ...
                      'psname','preprocess.ps','noiseth',p.noiseth,...
                      'loglevel',p.loglevel);
    spm12w_logger('msg','Slice noise check complete...', 'level',p.loglevel) 
    % Print then close slice noise figure
    print(F, 'preprocess.ps', '-dpsc2','-painters','-append','-noui')
end

% Step: Delete intermediate and unncessary files
if p.cleanup > 0
    % 0 = keep all files (useful for debugging)
    % 1 = delete all but the last stage of preprocessing (usually 's')
    % 2 = delete all but the last 2 stages of preprocessing (usually 'sw')
    % 3 = delete all but the last 2 stages of prerocessing and the raw epi
    deltmp = p.fmri(1:strfind(p.fmri,'epi')-1); % get only the prefixes.
    switch(p.cleanup)
        case 1
            keep_str = deltmp(1);            
            del_str = deltmp(2:end);
            d_msg = sprintf(['Cleanup (option:1), preserving ' ...
                            'the final stage of preprocessing: %s'], ...
                            keep_str); 
        case 2
            keep_str = deltmp(1:2);            
            del_str = deltmp(3:end);
            d_msg = sprintf(['Cleanup (option:2), preserving ' ...
                            'the final 2 stages of preprocessing: %s'], ...
                            keep_str);
        case 3
            keep_str = deltmp(1:2);
            del_str = deltmp(3:end);
            d_msg = sprintf(['Cleanup (option:3), preserving raw epi ' ...
                            'and the final 2 stages of preprocessing: %s'], ...
                            keep_str);
    end
    spm12w_logger('msg',d_msg,'level',p.loglevel) 
    % Iterate through del_str and create list of files to delete. Checking 
    % for whether we also delete the mat files and the raw epi files 
    % (i.e., p.matkill and p.cleanup ~= 3). 
    del_flist = {};
    for ses = 1:p.nses
        for i = 1:length(del_str)
            del_flist{end+1} = sprintf('%sepi_r%02d.nii',del_str(i:end),ses);
            if p.matkill
                if exist(fullfile(p.datadir,sprintf('%sepi_r%02d.mat',del_str(i:end),ses)),'file')
                    del_flist{end+1} = sprintf('%sepi_r%02d.mat',del_str(i:end),ses);
                end
                if exist(fullfile(p.datadir,sprintf('%sepi_r%02d_uw.mat',del_str(i:end),ses)),'file')
                    del_flist{end+1} = sprintf('%sepi_r%02d_uw.mat',del_str(i:end),ses);
                end
            end
        end
        if p.cleanup ~= 3
            del_flist{end+1} = sprintf('epi_r%02d.nii',ses);
        end     
    end
    % find the mean raw epi file whose exact name can vary depending on
    % preprocessing steps
    meanfile = dir(fullfile(p.datadir,'mean*epi_r01.nii'));
    if length(meanfile) == 1 && p.cleanup ~=3
        del_flist{end+1} = meanfile.name;
        if p.matkill
            [~,tmpname] = fileparts(meanfile.name);
            del_flist{end+1} = [tmpname,'_sn.mat'];
        end
    end
    % Delete every file in del_flist
    for del_f = del_flist
        spm12w_logger('msg',sprintf('[DEBUG] Removing file: %s', ...
                          del_f{1}),'level',p.loglevel)
        delete(fullfile(p.datadir,del_f{1}))
    end
    % If p.cleanupzip = 1, zip up the final stage files.
    if p.cleanupzip 
        zip_flist = {};
        for ses = 1:p.nses
            for i = 1:length(keep_str)
                zip_flist{end+1} = sprintf('%sepi_r%02d.nii',[keep_str(i:end),del_str],ses);
            end
            if p.cleanup == 3
                zip_flist{end+1} = sprintf('epi_r%02d.nii',ses);
            end   
        end
        % find the mean epi file whose exact name can vary depending on
        % preprocessing steps
        meanfile = dir(fullfile(p.datadir,'wmean*epi_r01.nii'));
        if length(meanfile) == 1
            zip_flist{end+1} = meanfile.name;
        end
        % Zip every file in del_flist and delete original
        for zip_f = zip_flist
            spm12w_logger('msg',sprintf('[DEBUG] Gzipping file: %s', ...
                          zip_f{1}),'level',p.loglevel)
            gzip(fullfile(p.datadir,zip_f{1}))
            delete(fullfile(p.datadir,zip_f{1}))
        end
    end        
    spm12w_logger('msg','Cleanup complete...', 'level',p.loglevel)     
end

% Close hidden figure (try because in some cases it might already be closed)
try
    F = spm_figure('FindWin','Graphics');
    close(F)
end

% Convert multipage ps file to pdf using spm12w_ps2pdf.m
spm12w_ps2pdf('ps_file',fullfile(p.datadir,'preprocess.ps'),...
              'pdf_file',fullfile(p.datadir,[p.sid,'_',p.prep_name,'.pdf']));

% Copy qa pdf to qa dir.
copyfile(fullfile(p.datadir,[p.sid,'_',p.prep_name,'.pdf']), ...
         fullfile(p.qadir,[p.sid,'_',p.prep_name,'.pdf']));

% Save parameter structure to mat file
save([p.prep_name,'.mat'],'p');
   
% Final words
msglist{1} = p.niceline;
msglist{2} = sprintf('Preprocessing complete on subject: %s',p.sid);
msglist{3} = sprintf('Log file   : %s', fullfile(p.datadir, p.preplog));
msglist{4} = sprintf('Figures    : %s', fullfile(p.datadir,[p.sid,'_',p.prep_name,'.pdf']));
msglist{5} = sprintf('Parameters : %s', fullfile(p.datadir,[p.prep_name,'.mat']));
if p.snr
    msglist{6} = sprintf('QA output  : %s', fullfile(p.qadir));
end
for msg = msglist
    spm12w_logger('msg',msg{1},'level',p.loglevel)
end

% Close log and return to studydir.
diary off; 
cd(p.study_dir)