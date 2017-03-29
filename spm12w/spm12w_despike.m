function spm12w_despike(varargin)
% spm12w_despike(niifile, cut, dilate)
%
% Input
% -----
% niifile  : Full path to nifti file for despiking. Can accept .gz files.  
%
% outfile  : Name of nifti file for output of despiking. 
%
% cut      : Cut range for despiking [c1,c2] (default=[2.5,4]).  
%
% c_order  : Curve fit order. Sets the number of sin and cosine terms to
%            fit. (Default is number of volumes / 30). 
%
% polort   : Order of polynomials to include (afni's 3dDespike uses 2)
%            0:Const, 1:Linear, 2:Quadratic, 3:Cubic (default=2).
%
% least_sq : Switch to enable least squares regression, otherwise proceeds
%            with afni 3dDespike default of L1 (Least Absolute Deviation) 
%            regression.
%
% dilate   : Amount to dilate an epi mask generated with spm's ??? function.
%            Dilating the mask will ensure despike is run on all epi voxels
%            but will maintain the speed boost of not running over the
%            entire volume space. 
%
% validate : If a voxel is provided (e.g., [22,22,22]) then L1 and OLS based
%            curve fits will be generated for that voxel along with a figure 
%            shwoing the raw data and the curve fits to that data. This is
%            to validate curve fitting and verify that everything is
%            working as expected.
%
% loglevel : Optional log level for spm12w_logger corresponding (default=1)
%
% Matlab implementation of afni's 3dDespike. 
%
% Examples:
%
% Despike a volume
%   >> spm12w_despike('niifile', 'epi_r01.nii', 'cut', [2.5,5], ...
%                     'c_order',5, 'polort',3, 'dilate', 4)
%
% Generate validation figures
%   >> spm12w_despike('niifile', 'epi_r01.nii', 'validate', [22,22,22])
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: November, 2014 | Updated: March, 2017
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('niifile','', 'outfile','despiked.nii', 'cut',[2.5,4], ...
              'polort',2, 'least_sq',0, 'c_order','nvols/30', ...
              'dilate',4, 'validate', '', 'loglevel', 1);
args = spm12w_args('nargs', 2, 'defaults', args_defaults, 'arguments', varargin);

% Check if validate is empty. 
% If empty, proceed as usual. 
% Otherwise generate curve fitting validation figure. 
if isempty(args.validate)
    niidata = spm12w_readnii('niifiles', args.niifile, 'loglevel', args.loglevel);
    % Generate despike design matrix (xmat) with internal function
    xmat = make_xmat(length(niidata.hdrs), args.c_order, args.polort); 
    % Create empty data matrix
    ds_data = zeros(size(niidata.data));
    % loop over voxels (can we do this more effeicently? slicewise?)
    for i_3 = 1:size(niidata.data,3)        
        spm12w_logger('msg',sprintf('Despiking: slice %d of %d', ...
              i_3, size(niidata.data,3)), 'level',args.loglevel)
        for i_2 = 1:size(niidata.data,2)
            for i_1 = 1:size(niidata.data,1)
                y = squeeze(niidata.data(i_1,i_2,i_3,:));
                % Compute regressions 
                if args.least_sq == 1
                    b = xmat \ y;
                else %l1 regression
                    b = spm12w_stats('stat', 'l1', 'y', y, 'x', xmat);
                end
                % Identify spikes by looking at deviation from curve for
                % each timepoint.
                yhat = xmat * b;      % Get predicted values
                % ymad = median-mad of abs value of deviation from curve.
                ymad = spm12w_stats('stat', 'mad_med', ...
                                    'y', abs(y-yhat)); 
                fsig = sqrt(pi/2)*ymad;   % stdev of residuals
                yspike = (y-yhat)/fsig;   % identify spikes
                % replace spikes with formula: s' = c1+(c2-c1)*tanh((s-c1)/(c2-c1))
                for i_sp = 1:length(yspike)
                    if yspike(i_sp) > args.cut(1)
                        %fprintf('%.4f',yspike(i_sp))
                        yspike(i_sp) = args.cut(1)+(args.cut(2)-args.cut(1))*tanh((yspike(i_sp)-args.cut(1))/(args.cut(2)-args.cut(1)));
                        %fprintf('...%.4f\n',yspike(i_sp))
                    elseif yspike(i_sp) < -args.cut(1)
                        %fprintf('%.4f',yspike(i_sp))
                        yspike(i_sp) = -(args.cut(1)+(args.cut(2)-args.cut(1))*tanh((abs(yspike(i_sp))-args.cut(1))/(args.cut(2)-args.cut(1))));
                        %fprintf('...%.4f\n',yspike(i_sp))
                    end
                end
                ds_data(i_1,i_2,i_3,:) = yspike;
            end
        end
    end
    % Write despiked data to output filename
    spm12w_writenii('niifile',args.outfile, ...
                    'hdr',niidata.hdrs,'data',ds_data);                     
else
    % Validation of curve fitting technique against afni's 3dDespike curve
    % --------------------------------------------------------------------
    % This appears accurate in that we're fitting an interesting curve to the
    % data. The formula I've used appears to be identical to that
    % in the afni 3dDespike docs. And the L1 regression solution should be 
    % approximatly the same. True validation will require comparing voxels from
    % 3dDespike and this method.

    % Extract timeseries at specifed validation voxel from the supplied niifile
    spm12w_logger('msg', sprintf(['Running despike validation on voxel: ', ...
                  '%d,%d,%d'], args.validate), 'level',args.loglevel)
    y = spm12w_readnii('niifile',args.niifile, 'voxels',args.validate, ...
                       'vox_only', 1);
    y = y.voxdata;
    xmat = make_xmat(length(y), args.c_order, args.polort); %generate xmat with defaults
    % L1 Regression
    bl1 = spm12w_stats('stat', 'l1', 'y', y, 'x', xmat);
    yhat_l1 = xmat * bl1; % Get predicted values
    % Least squares regression 
    bls = xmat \ y; % equivalent inv(xmat'*xmat)*xmat'*y; 
    yhat_ls = xmat * bls; % Get predicted values

    % Curve fitting toolbox method
        % In practice this has always been identical to the least squares fits.
        % ft = fittype( 'aa+bb*x+cc*x^2+dd*sin(2*pi*x*1/120)+ee*cos(2*pi*x*1/120)+ff*sin(2*pi*x*2/120)+gg*cos(2*pi*x*2/120)+hh*sin(2*pi*x*3/120)+ii*cos(2*pi*x*3/120)+jj*sin(2*pi*x*4/120)+kk*cos(2*pi*x*4/120)', 'independent', 'x', 'dependent', 'y' );
        % opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        % opts.Display = 'Off';
        % opts.StartPoint = [0.259902669487639 0.684443374772421 0.260210359730971 0.129067975009114 0.101872152409411 0.975009916772407 0.251162172997024 0.0994406427539448 0.63685208690899 0.224963657822537 0.81067965939345];
        % [fitresult, gof] = fit( x, y, ft, opts );
        % bcf = coeffvalues(fitresult)';
        % yhat_cf = xmat * bcf;

    % Plot figure of L1 and LS methods.
    ts_label = sprintf('Timeseries data (%d,%d,%d)',args.validate);
    figure( 'Name', 'Validation of spm12w_despike curve fitting' );
    X = 1:length(y);
    plot(X,y,'m:o',...
        X,yhat_l1,'r-',...
        X,yhat_ls,'k-')
    axis square
    grid on
    legend(ts_label,'L-1 Regression','Least Squares Regression')
end

function xmat = make_xmat(nvols, c_order, polort)
% Build a design matrix of predictors consisting of:
% Constant, Polynomial Trends, sin and cos trends.
% Note: we can't use linspace(-1,1) for polys since it makes it difficult 
%       to plot curves and get y^hat in same units as data. This method is 
%       equivalent (as far as I can tell) to the afni C source code. 
    if isinteger(c_order)
        c = c_order;               % Curve fit order (see afni's 3dDespike);
    else
        c = round(eval(c_order));
    end
    polys = 0:polort;              % Full set of polys
    n_pred = 2*c+polort+1;         % Number of predictors
    x = (1:nvols)';                % x-values are n_vols.
    xmat = zeros(nvols, n_pred);   % Zero design matrix 

    for i = 1:length(polys)
        xmat(:,i) = x.^polys(i); % Add polys
    end
    % Create predictors for sin and cos up to corder. 
    i_pred = i+1; % Set current predictor to i+1 (i.e., the last poly above)
    for i_c = 1:c
        xmat(:,i_pred) = sin(((1:nvols)'*2*pi*i_c)/nvols);
        i_pred = i_pred + 1;
        xmat(:,i_pred) = cos(((1:nvols)'*2*pi*i_c)/nvols);
        i_pred = i_pred + 1;
    end
return
