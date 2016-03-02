function spm12w_snr(varargin)
% spm12w_snr(niifile, outdir, sid, psname, detrend, loglevel)
%
% Input
% -----
% niifiles : A cell array of full paths to nifti files(s) for tSNR
%            calcualtion.
%
% outdir   : Path to directory where tSNR niftifiles and QA figures will be
%            written. If none is specified, results will be written to current
%            directory. 
%
% sid      : Subject ID prefix for output nifti files (e.g., snr, sd, mean)
%            If no sid is specified, results will be written without sid prefix. 
%
% psname   : Name of .ps file containing figures displaying SNR plots.
%            Default is preprocess.ps
%
% detrend  : Flag for detrending data prior to tSNR calculation. tSNR may
%            potentially be inflated by trends in the data that are usually
%            removed during GLM estimation. This switch will do detrending 
%            prior to tSNR calcualtion. This adds considerable processing time. 
%            (default=0) (AND ISN'T ACTUAL IMPLEMENTED!)
%
% scalefig : Flag for scaling data to fit in the display range for showing
%            figures. This isn't perfect but helps with visualization. The
%            actual nifti files that are written out are not scaled. This
%            addresses the issue of different scanners outputting values not
%            in the range that this script exepcts. (default=1). 
%
% loglevel : Optional log level for spm12w_logger corresponding (default=1)
%
% Calculates mean image, sigma and temporal SNR (i.e., tSNR = mean/sigma).
% Creates figures and outputs qa files to the directory specified in output
% directory. 
%
% Examples:
%
% Calculate snr on two runs of data
%   >> spm12w_snr('niifiles', {'./epi_r01.nii', './epi_r02.nii'})
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: November, 2014 | Updated: March, 2016
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('niifiles','','outdir',pwd,'sid','',...
                       'psname','preprocess.ps','detrend',0,...
                       'scalefig',1,'loglevel', 1);
args = spm12w_args('nargs',2, 'defaults', args_defaults, 'arguments', varargin);

% Make outdir if it doesn't exist.
if ~isdir(args.outdir)
    mkdir(args.outdir)
end

% Setup figure (if it exists, clear it, otherwise make it).
F = spm_figure('FindWin','Graphics');
if isempty(F)
    F = spm_figure('FindWin','Graphics', 'spm12w preprocessing', 'off'); 
else
    spm_figure('Clear',F);
end

% Loop through files (i.e., sessions) making tsnr calculations and plots
Vnii = spm_vol(args.niifiles);
currentplot = 1; % sublot for figure
for ses = 1:size(Vnii,2)
    % Check if only one session/niifile (in which case Vnii is struct not cell)
    if iscell(Vnii)
        files  = Vnii{ses};
        spm12w_logger('msg',sprintf(['[DEBUG] Calculating tSNR for ' ...
                      'session %d of %d'], ses, size(Vnii,2)), ...
                      'level',args.loglevel)      
    else 
        files = Vnii;
    end
    % Get scaling factor
    scalf  = files(1).private.dat.scl_slope;
    if scalf > 0
        spm12w_logger('msg',sprintf(['[DEBUG] Data will be scaled by ' ...
                      '%.2f prior to SNR calculation'], scalf), ...
                      'level',args.loglevel) 
    end
    % Setup tsnr variables.
    data    = spm_read_vols(files(1));
    avg     = zeros(size(data));
    delta   = zeros(size(data));
    datavar = zeros(size(data));
    nvols   = size(files,1);
    % Incremental algorithm for rolling mean and variance
    % (saves memory by not preloading all volumes). 
    for i = 1:nvols          
        data = spm_read_vols(files(i))*scalf; %scale data prior to calc
        delta = data - avg;
        avg = avg+delta/i;
        datavar = datavar + delta.*(data-avg);
    end
    sd = sqrt(datavar/(nvols-1));
    snr = avg./sd;
    % Clean up snr varaible (remove nan and ridiculously high values (due
    % to low variance). 
    snr(isnan(snr)) = 0;
    snr(snr>5000)   = 0; 
    % Output avg, sd and snr volumes
    spm12w_logger('msg','[DEBUG] Writing out volumes...', 'level',args.loglevel) 
    vol_out       = files(1);
    vol_out.pinfo = [1 0]'; %set scaling to 1 (i.e., we don't orig data scaling)
    vol_out.fname = fullfile(args.outdir,sprintf('avg_%s_r%02d.nii', args.sid, ses));
    spm_write_vol(vol_out, avg);
    vol_out.fname = fullfile(args.outdir,sprintf('sd_%s_r%02d.nii', args.sid, ses));
    spm_write_vol(vol_out, sd);
    vol_out.fname = fullfile(args.outdir,sprintf('snr_%s_r%02d.nii', args.sid, ses));
    spm_write_vol(vol_out, snr);
    % Make figures
    % Normalize average and SD. SNR should not be averaged as it is not scale
    % dependent. We'll use the same scale for the SD as the average. Seems to
    % work better when we do that.
    if args.scalefig
        spm12w_logger('msg','[DEBUG] Scaling for figure display purposes.', ...
                      'level',args.loglevel) 
        scaleavg = max(avg(:))/900;
        avg = avg/scaleavg;
        sd = sd/scaleavg;        
    end
    % Slices (assumes normalized data, these are decent spots)
    slice{1} = squeeze(avg(:,:,20));   stitles{1} = 'Average 1';   sthresh{1} = [10,900];
    slice{2} = squeeze(avg(:,:,24));   stitles{2} = 'Average 2';   sthresh{2} = [10,900];
    slice{3} = squeeze(sd(:,:,20));    stitles{3} = 'Stdev 1';     sthresh{3} = [2,20];
    slice{4} = squeeze(sd(:,:,24));    stitles{4} = 'Stdev 2';     sthresh{4} = [2,20];
    slice{5} = squeeze(snr(:,:,20));   stitles{5} = 'SNR 1';       sthresh{5} = [10,350];
    slice{6} = squeeze(snr(:,:,24));   stitles{6} = 'SNR 1';       sthresh{6} = [10,350];
    slice{7} = squeeze(snr(26,:,:));   stitles{7} = 'SNR Sagital'; sthresh{7} = [10,350];
    slice{8} = squeeze(snr(:,32,:));   stitles{8} = 'SNR Coronal'; sthresh{8} = [10,350];
    % Plot defaults
    s_xpos = [0.05,0.26,0.54,0.75,0.05,0.26,0.54,0.75];
    s_tpos = [27,27,27,27,27,27,33,27];
    colormap hot
    % Set Currentplot defaults (1 = TOP)
    if currentplot == 1 %Top plot
        t_ypos = 0.75;
        splot = 1;
        s_ypos = [0.75,0.75,0.75,0.75,0.55,0.55,0.55,0.55];
    else                %Bottom plot
        t_ypos = 0.27;
        splot = 9;
        s_ypos = [0.27,0.27,0.27,0.27,0.07,0.07,0.07,0.07];    
    end
    % Plots: 8 per halfpage. Subplot is 4,4, and 1 to 16.
    for iplot = 1:8
        subplot(4,4,splot)
        set(gca,'position',[s_xpos(iplot),s_ypos(iplot),0.2,0.2])
        imagesc(flipud(slice{iplot}'),sthresh{iplot})
        axis equal 
        axis off
        title(stitles{iplot},'fontweight','bold','position',[s_tpos(iplot),0.5])
        splot = splot + 1;
    end
    % Titles
    titlestr = ['SNR Subject: ',args.sid, ' Run: ',num2str(ses)];
    titleax = axes('Position',[0.1 t_ypos 0.8 0.2],'Parent',F,'Visible','off');
    set(get(titleax,'Title'),'String',titlestr,'FontSize',16,'FontWeight','Bold','Visible','on');
    % Print Check
    if ses == size(Vnii,2)
       print(F, args.psname, '-dpsc2','-painters','-append','-noui')       
       break;
    end            
    if currentplot == 1
        currentplot = 2;
    else
        print(F, args.psname, '-dpsc2','-painters','-append','-noui')      
        spm_figure('clear',F);  
        currentplot = 1;
    end     
end