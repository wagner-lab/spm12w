function spm12w_slicenoise(varargin)
% spm12w_slicenoise(niifile, radata, psname, loglevel)
%
% Input
% -----
% niifiles : A cell array of full paths to nifti files(s) for tSNR
%            calcualtion.
%
% radata   : A cell array of realignment data (i.e., load realignment 
%            parameters into an array or use p.ra if called from
%            spm12w_preprocess.m)
%
% mask     : Full path to a mask file to constrain slice calcualtions
%            (usually the default mask contained in p.mask)
%
% psname   : Name of .ps file containing figures displaying SNR plots.
%            Default is preprocess.ps
%
% noiseth  : Noise thresholds for displaying figures. Might need to be tweaked 
%            according to scanner. (default=[0.015,0.03]).
%
% loglevel : Optional log level for spm12w_logger corresponding (default=1)
%
% Calculates slice-wise differences (current - previous slic) and creates a 
% figure showing changes in overall slice mean. Borrowed from slices_analyse
% by Antonia Hamilton circa 2006. Data is mean centered slice-wise to help
% account for scanner / data scaling differences, not sure if this is ideal
% but works well enough.
%
% Examples:
%
% Calculate slice noise on two runs of data
%   >> spm12w_slicenoise('niifiles', {'./epi_r01.nii', './epi_r02.nii'},
%                        'radata', p.ra, 'mask', p.mask, 'noiseth',[0.015,0.03])
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: November, 2014 | Updated: May, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('niifiles','','radata','','mask','',...
                       'noiseth',[5,15],'psname','preprocess.ps',...
                       'loglevel', 1);
args = spm12w_args('nargs',4, 'defaults', args_defaults, 'arguments', varargin);

% Setup figure (if it exists, clear it, otherwise make it).
F = spm_figure('FindWin','Graphics');
if isempty(F)
    F = spm_figure('GetWin','Graphics', 'spm12w preprocessing', 'off'); 
else
    spm_figure('Clear',F);
end
  
% Load runs and concat into single run. 
V      = spm_vol(args.niifiles);
V      = cat(1,V{:}); 
ra     = cat(1,args.radata{:});
nslice = V(1).dim(3);
nscan  = length(V);
scalf = zeros(nscan,1);
for i = 1:nscan
    scalf(i) = V(i).private.dat.scl_slope;
end
if scalf(1) > 0
    spm12w_logger('msg',['[DEBUG] Data will be scaled prior to slice ' ...
              'noise calculation'],'level',args.loglevel)
end
% This assumes data is normalized.
mask_name = args.mask;
M         = spm_vol(mask_name);
M.dat     = spm_read_vols(M);
% work out good neighbours
neigh              = repmat(-5:2:5,nscan,1)+repmat(1:nscan,6,1)';
neigh(neigh<1)     = NaN;
neigh(neigh>nscan) = NaN;
for kk=1:nscan      
    % calculate neighbours to pull from for noise calc
    df      = abs(ra-repmat(ra(kk,:),nscan,1));
    trans   = sum(df(:,1:3),2);
    rot     = sum(df(:,4:6),2);
    good    = find(trans<1 & (rot*180/pi)<1);
    overlap = intersect(good,neigh(kk,:)); 
    if(length(overlap)<3)      % really bad scans
        neigh(neigh==kk) = NaN;     
        neigh(kk,:)      = NaN;
    elseif(length(overlap<6))  % moderate bad scans
        missing          = setdiff(neigh(kk,:),good);
        for i=1:length(missing)
            mind           = find(neigh(kk,:)==missing(i));
            neigh(kk,mind) = NaN;
        end
    end
end
% Calculate slice noise
spm12w_logger('msg',sprintf(['[DEBUG] Calculating slice noise for %d ', ...
              'volumes'], nscan),'level',args.loglevel)
for kk=1:nscan  % for every input file
    gv(kk) = spm_global(V(kk));
    [dat1,loc] = spm_read_vols(V(kk),1);  % read with zero masking 
    dat1 = dat1/gv(kk);         % Apply scaling
    nn = neigh(kk,:);
    for jj=1:length(nn)    % read some neighbours
        if(isfinite(nn(jj)))  % for good neighbours
            jind = nn(jj);
            gv2(jj) = spm_global(V(jind));
            [dat2,loc] = spm_read_vols(V(jind),1);  %% read with zero masking
            dat2 = dat2/gv2(jj);                %% Apply scaling
            for i=1:nslice
                slice1 = squeeze(dat1(:,:,i));
                slice2 = squeeze(dat2(:,:,i));
                msk = squeeze(M.dat(:,:,i));
                df = (slice1-slice2).*(msk>0.5);   %mask the data      
                scan_noise(jj,i) = nanmean(abs(df(msk>0.5)));
            end
        else  % for bad neighbours don't calc slicenoise
            scan_noise(jj,1:nslice) = NaN;
        end
    end
    % Average the mean of the difference between current slice
    % and all good slice neighbors
    noise(:,kk) = nanmean(scan_noise)'; 
end
% Make slices figure (now 200% more better!!!)
th  = args.noiseth(1);  %default was 20 from slices_defaults.m
wth = args.noiseth(2);  %default was 30 lowering defaults since new coil -DDW Feb/12 
noise = noise*400; % Arbitrary hack to get thresholds in a good spot. 
colormap('default');
subplot(3,1,1)
imagesc(noise,[0,40])
ttl = ['Filetype: ',strtok(spm_str_manip(args.niifiles{1},'t'),'_')];
title(ttl);
subplot(3,1,2)
[n,b] = hist(noise(:),[0:40]);
bar(b,n);
set(gca,'xlim',[0,40])
title('distribution of slice noise')
hold on
plot([th,th],[0,max(n)],'r-')
h=text(th,max(n),[num2str(100*sum(noise(:)>th)./numel(noise),3),...
              '% of slices are over ',num2str(th)]);
set(h,'Color',[1 0 0])
plot([wth,wth],[0,max(n)/2],'g-')
h=text(wth,max(n)/2,[num2str(100*sum(noise(:)>wth)./numel(noise),3),...
                  '% of slices are over ',num2str(wth)]);
set(h,'Color',[0 1 0])
subplot(6,1,5)
plot(ra(:,1:3))
title('translation (mm)')
set(gca,'XTick',[])
subplot(6,1,6)
plot(ra(:,4:6)*180/pi)
title('rotation (deg)')
% Title
titleax = axes('Position',[0.12 0.75 0.8 0.2],'Parent',F,'Visible','off');
set(get(titleax,'Title'),'String','Slice Noise Analysis','FontSize',16,'FontWeight','Bold');
% Print
spm12w_logger('msg',['[DEBUG] Printing slice noise figure'],'level',args.loglevel)
print(F, args.psname, '-dpsc2','-painters','-append','-noui')       