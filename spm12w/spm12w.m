function spm12w(varargin)
% spm12w(stage, sids, para_file, parloop, parcores)
%
% Inputs
% ------
% stage: String for the specific analysis stage desired. Supplying a stage
%        will override the gui and go directly to performing the analysis
%        stage. If stage requires a p & glm parameters file you must provide
%        both as cells in the para_file inputs. 
%        Stages are: prep, prep_glm, prep_glm_con, glm, glm_con, rfx, roi. 
%
% sids: Cell array of subject ids. If left unspecified, a dialog box will
%       appear asking the user to select subjects. If stage is empty, then
%       specifying sids will not override the gui. You may also specify
%       the path to a .mat file containing a saved sidlist cell array 
%       (e.g. sidlist = {'s01','s04','s12'}; save('sidlist.mat','sidlist'))
%
% para_file: If using spm12w in non-gui mode a cell array of para_files is
%            requires (e.g., para_file = {'./scripts/username/p_tutorial.m'})
%
% parloop: Optional input (0 or 1) to enable parallel processing, useful when
%          overriding the gui. Default is no parallel processing (default = 0). 
%
% parcores: Optional input to specify the size of the parallel processing pool.
%           Can be set in gui, but useful to specify here if overriding the gui.
%           Default value is 4 cpu cores but can range from 2 to 64. Caution
%           should be used when enabling this on a shared server environment.
%           
% Examples:
% Please note that spm12w should always be run from a study's top-level 
% directory.
%
% To select analyses via the spm12w gui, type spm12w with no argument:
%    
%       >> spm12w
%
% To run spm12w preprocessing on sids s01 and s02 overriding the gui:
%  
%       >> spm12w('stage','prep','sids',{'s01','s02'},...
%                 'para_file',{'./scripts/username/p_tutorial.m'})
%
% To run spm12w preprocessing on sids s01 and s02 w/o overriding the gui:
%  
%       >> spm12w('sids',{'s01','s02'})
%
% To run spm12w preprocessing on a saved list of sids overriding the gui:
%  
%       >> spm12w('stage','prep','sids','./auxil/sidlist.mat',...
%                 'para_file',{'./scripts/username/p_tutorial.m'})
%
% To run spm12w preprocessing and GLM on a saved list of sids and overriding gui
%
%       >> spm12w('stage','prep','sids','./auxil/sidlist.mat',...
%                 'para_file',{'./scripts/username/p_tutorial.m',...
%                 './scripts/username/glm_tutorial.m'})
%
% To run spm12w preprocessing and GLM on a saved list of sids and overriding gui
% and providing parallele processing argumetns (i.e, the works!). 
%
%       >> spm12w('stage','prep','sids','./auxil/sidlist.mat',...
%                 'para_file',{'./scripts/username/p_tutorial.m',...
%                 './scripts/username/glm_tutorial.m'}, ...
%                 'parloop',1,'parcores',24)
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: May, 2015 | Updated: March, 2016
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('stage','', 'sids','','para_file','',...
                       'parloop',0, 'parcores', 4);
args = spm12w_args('nargs',0, 'defaults', args_defaults, 'arguments', varargin);

% Clear functions is necessary to force update of .m files that have been
% changed. Ocassionally Matlab fails to notice an updated parameters (.m)
% file and runs a prior version cached in memory. Seems to affect linux 
% primarily but no harm running on pc/mac
clear functions

% assign cwd for oncleanup
cwd=pwd;

% assign global gui data structure for passing data between functions. Later
% consider not using a global variable, but for now it's the simplest way
% to pass data between callbacks and functions.
global gui
% Set defaults
gui.niceline = [repmat('=',1,80)];
gui.name = 'spm12w-gui';
gui.imgfile = which('spm12w_title.jpg');
% Define color pallet (push this to defaults later)
gui.p_bcknd = [206,198,189] ./255;
gui.p_panel = [80,81,79] ./255;
gui.p_ui = [.64,.64,.64];
gui.p_uitxt = [255,204,102] ./255;
gui.p_border = [20,0,0] ./255;
gui.p_title1 = [255,224,102] ./255;
gui.p_title2 = [242,125,122] ./255;
gui.p_title3 = [222,110,75] ./255;
% Assign variables for parallel processing
gui.parloop = args.parloop;
gui.parcores = args.parcores;

% Check for user using strings instead of cells.
if ischar(args.sids) && ~isempty(args.sids)
    args.sids = cellstr(args.sids);
end

% Check of user supplying a string para_file instead of a cell array
if ischar(args.para_file) && ~isempty(args.para_file)
    args.para_file = cellstr(args.para_file);
end

% Check for user supplied sidlist mat file
if length(args.sids) == 1 && any(strfind(args.sids{1},'.mat'))
    spm12w_logger('msg',sprintf('[DEBUG] Loading sidlist .mat file: %s', ...
                                 args.sids{1}));
    args.sids = load(args.sids{1});
    args.sids = args.sids.sids; %fix for matlab loading into struct array
end

% If user submitted an analysis stage then do not invoke gui.
if isempty(args.stage)
    gui.nogui = 0;
    ressurect_gui(args.sids, args.para_file) % Arise you shinning GUI!
else
    gui.nogui = 1;
    spm12w_logger('msg',gui.niceline)
    spm12w_logger('msg',sprintf('Overriding gui and running stage: %s', ...
                             args.stage))
    if gui.parloop == 1
        spm12w_logger('msg',sprintf(['[DEBUG] Parallel processing is ',...
                      'enabled using %d cores:'], gui.parcores))
    end
    stager([],[],args.stage, args.sids, args.para_file) 
end
  
% on cleanup return to root.
croot = onCleanup(@()cd(cwd));

function ressurect_gui(sids, para_file)
   global gui   
   close all
   % Define gui positions
   gui.pos_gui = [50,400,552,658]; %Position to place gui figure (X,Y,W,H)
   % Generate gui
   gui.spm12wfig = figure('Name',gui.name, ...
                     'MenuBar','none', ...
                     'Toolbar','none', ...
                     'NumberTitle','off', ...
                     'Resize', 'off', ...
                     'Visible', 'off', ...
                     'Position',gui.pos_gui, ...
                     'color', gui.p_bcknd); 
   % Set gui image
   titleimage = axes('Units','pixels','position',[2,581,550,107]); 
   image(imread(gui.imgfile)); 
   axis('off','image') 
   % Define Panels
   ypos = 524; %Sets starting Y position. All other is relative to this
   % Main panel
   p0 = uipanel('BackgroundColor',(gui.p_panel),'Units','pixels',...
                'Position',[52,6,450,600],'BorderWidth',2,...
                'BorderType','line','HighlightColor',gui.p_border); 
   % Preprocessing Panel
   p1 = uipanel('BackgroundColor',(gui.p_panel),'Units','pixels',...
                'Position',[60,ypos,435,78],'title','Preprocessing',...
                'HighlightColor','black','BorderType','line',...
                'ForegroundColor',gui.p_title1); 
   % GLM / RFX Panel /ROI 
   p2 = uipanel('BackgroundColor',(gui.p_panel),'Units','pixels',...
                'Position',[60,(ypos-82),435,78],'title',...
                '1st & 2nd Level Analysis','HighlightColor','black',...
                'BorderType','line','ForegroundColor',gui.p_title1); 
   % ROI & Conjunction
   p3 = uipanel('BackgroundColor',gui.p_panel,'Units','pixels',...
                'Position',[60, (ypos-144),435,58],'title',...
                'ROI & Conjunction Analysis','HighlightColor','black',...
                'BorderType','line','ForegroundColor',gui.p_title1); 
   % Segment and DARTEL Normalization
   p4 = uipanel('BackgroundColor',gui.p_panel,'Units','pixels',...
                'Position',[60, (ypos - 226),435,78],'title',...
                'Segment & DARTEL Normalization','HighlightColor','black',...
                'BorderType','line','ForegroundColor',gui.p_title1);
   % VOI & PPI Analysis
   p5 = uipanel('BackgroundColor',gui.p_panel,'Units','pixels',...
                'Position',[60, (ypos-288),435,58],'title',...
                'VOI Extraction & PPI Regressor Maker', ...
                'HighlightColor','black','BorderType','line', ...
                'ForegroundColor',gui.p_title1);            
   % Design search & build
   p6 = uipanel('BackgroundColor',gui.p_panel,'Units','pixels',...
                'Position',[60, (ypos - 370),435,78],'title',...
                'Experiment Design Helper','HighlightColor','black',...
                'BorderType','line','ForegroundColor',gui.p_title2);         
   % Utilities
   p7 = uipanel('BackgroundColor',gui.p_panel,'Units','pixels',...
                'Position',[60, (ypos - 452),435,78],'title',...
                'Utilities','HighlightColor','black',...
                'BorderType','line','ForegroundColor',gui.p_title2);
   % Parfor switch
   p8 = uipanel('BackgroundColor',gui.p_panel,'Units','pixels',...
                'Position',[60, (ypos-514),435,58],'title',...
                'Matlab Parallel Computing Toolbox',...
                'HighlightColor','black','BorderType','line',...
                'ForegroundColor',gui.p_title3);    
   % Define stage text 
   stage={'Preprocess','prep';
          'Preprocess & GLM','prep_glm';
          'Preprocess & GLM & Contrasts','prep_glm_con';
          'Estimate GLM','glm';
          'Compute Contrasts','con';
          'Estimate GLM & Con','glm_con';           
          'Random Effects','rfx';
          'ROI Analysis','roi';
          'WIP-Conjunction','conj';
          'WIP-Segment Anatomy','seg8';
          'WIP-Template Maker','tempflow';
          'WIP-Normalize Anat','nanat';
          'WIP-Normalize EPI','nfunc';
          'WIP-VOI Extraction','voi';
          'WIP-PPI Reg Maker','ppi';             
          'Design Search','dsearch';
          'Design Check','dcheck';
          'Design Build','dbuild'  
          'ART Outlier Detect','art';
          'Onset Maker','onsmk';
          'Prep Study','prepare';
          'WIP-Cluster Sims','alphasim'};
    radio={'[ON] Parallel';'[OFF] Serial';'CPU Cores:'};
    %Panel 01-Preprocessing
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{1,1},'Units', 'pixels','Position',[44,34,160,26],'Parent',p1,'Callback',{@stager, 'prep',sids, para_file});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{2,1},'Units', 'pixels','Position',[229,34,160,26],'Parent',p1,'Callback',{@stager,'prep_glm',sids, para_file});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{3,1},'Units', 'pixels','Position',[84,5,260,26],'Parent',p1,'Callback',{@stager, 'prep_glm_con',sids, para_file});
    %Panel 02-1st & 2nd Level Analysis
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{4,1},'Units', 'pixels','Position',[44,34,160,26],'Parent',p2, 'Callback',{@stager, 'glm',sids, para_file});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{5,1},'Units', 'pixels','Position',[229,34,160,26],'Parent',p2, 'Callback',{@stager, 'con',sids, para_file});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{6,1},'Units', 'pixels','Position',[44,5,160,26],'Parent',p2, 'Callback',{@stager, 'glm_con',sids, para_file});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{7,1},'Units', 'pixels','Position',[229,5,160,26],'Parent',p2,'Callback',{@stager, 'rfx',sids, para_file});    
    %Panel 03-ROI & Conjunction Analysis
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{8,1},'Units', 'pixels','Position',[44,10,160,26],'Parent',p3, 'Callback',{@stager, 'roi',sids, para_file});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{9,1},'Units', 'pixels','Position',[229,10,160,26],'Parent',p3,'Callback',{@stager, 'conj',sids});
    %Panel 04-Segment and DARTEL
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{10,1},'Units', 'pixels','Position',[44,34,160,26],'Parent',p4, 'Callback',{@stager, 'glm',sids, para_file});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{11,1},'Units', 'pixels','Position',[229,34,160,26],'Parent',p4, 'Callback',{@stager, 'con',sids, para_file});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{12,1},'Units', 'pixels','Position',[44,5,160,26],'Parent',p4, 'Callback',{@stager, 'glm_con',sids, para_file});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{13,1},'Units', 'pixels','Position',[229,5,160,26],'Parent',p4,'Callback',{@stager, 'rfx',sids, para_file});    
    %Panel 05-VOI Extraction & PPI Regressor Maker
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{14,1},'Units', 'pixels','Position',[44,10,160,26],'Parent',p5, 'Callback',{@stager, 'voi',sids});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{15,1},'Units', 'pixels','Position',[229,10,160,26],'Parent',p5,'Callback',{@stager, 'ppi',sids});
    %Panel 06-Experiment Design Helper 
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{16,1},'Units', 'pixels','Position',[44,34,160,26],'Parent',p6, 'Callback',{@stager, 'dsearch','',para_file});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{17,1},'Units', 'pixels','Position',[229,34,160,26],'Parent',p6, 'Callback',{@stager, 'dcheck','',para_file});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{18,1},'Units', 'pixels','Position',[134,5,160,26],'Parent',p6, 'Callback',{@stager, 'dbuild','',para_file});   
    %Panel 07-Utilities
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{19,1},'Units', 'pixels','Position',[44,34,160,26],'Parent',p7, 'Callback',{@stager, 'art',sids});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{20,1},'Units', 'pixels','Position',[229,34,160,26],'Parent',p7, 'Callback',{@stager, 'onsmk',sids});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{21,1},'Units', 'pixels','Position',[44,5,160,26],'Parent',p7, 'Callback',{@stager, 'prepare',sids});
    uicontrol('Style','PushButton','BackgroundColor', gui.p_ui,'HorizontalAlignment','left','String',stage{22,1},'Units', 'pixels','Position',[229,5,160,26],'Parent',p7,'Callback',{@stager, 'alphasim',sids});    
    %Set radio buttons for parfor
    par_bg = uibuttongroup(p8,'Units','pixels','Position',[10,8,420,30],'BorderType','line','HighlightColor',gui.p_border);
    uicontrol(par_bg,'Style','text','String',radio{3},'Position',[230,5,100,16],'tag','parcoretext');
    handles.cores = uicontrol(par_bg,'Style','edit','String',gui.parcores,'Position',[322,4,40,20],'tag','parcores','Callback',{@radio_clbk});
    uicontrol(par_bg,'Style','radiobutton','String',radio{1},'Position',[6,6,110,16],'Value',0,'tag','paron','Callback',{@radio_clbk, handles});
    uicontrol(par_bg,'Style','radiobutton','String',radio{2},'Position',[120,6,110,16],'Value',1,'tag','paroff','Callback',{@radio_clbk, handles}); 
    setcores = uicontrol(par_bg,'Style','PushButton','Tag','setcores','HorizontalAlignment','left','String','Set','Units','pixels','Position',[366,2,50,26],'Callback',{@radio_clbk, handles});    
    
    %Activate figure
    set(gui.spm12wfig, 'Visible', 'on') % Render the invisible visible.     
return

function radio_clbk(h, eventdata, handles)
    global gui
    switch get(h,'tag')
        case 'paron'
            gui.parloop = 1;
            spm12w_logger('msg',sprintf(['[DEBUG] Parallel ',...
                          'processing is [ON]']))                      
        case 'paroff'
            gui.parloop = 0;
            spm12w_logger('msg',sprintf(['[DEBUG] Parallel ',...
                          'processing is [OFF]']))       
        case 'parcores'
            gui.parcores = round(str2num(get(h,'String')));
        case 'setcores'
            gui.parcores = round(str2num(get(handles.cores,'string')));
            spm12w_logger('msg',sprintf(['[DEBUG] Parallel ',...
                          'processing set to %d cores...'],gui.parcores))          
    end    
return

% Implement the requested analysis stage
function stager(h, eventdata, rstage, sids, para_file)
    global gui  
        %Close spm12wfig as no longer needed 
    if ~gui.nogui
        close(gui.spm12wfig);
    end

    % Load sids and para_files for stages that require them
    sid_stages = {'prep','prep_glm','prep_glm_con','glm','glm_con','con',...
                 'rfx','roi'};
    p_stages   = {'prep','prep_glm','prep_glm_con'};
    glm_stages = {'prep_glm','prep_glm_con','glm','glm_con'};
    con_stages = {'prep_glm_con','glm_con','con'};        
    rfx_stages = {'rfx'};
    roi_stages = {'roi'};
    design_stages = {'dsearch','dcheck','dbuild'};
    
    % Load sids for sid_stages
    if ismember(rstage,sid_stages)
        if ismember(rstage,sid_stages(1:3))
            sids = spm12w_getsid(fullfile(pwd,'raw'));
        elseif ismember(rstage,sid_stages(4:end))
            sids = spm12w_getsid();
        end
        if isempty(sids)
            error('No subjects selected...');
        end  
    end
    
    % Load para files for p_stages, glm_stages, roi_stages or design_stages
    if isempty(para_file)
        if ismember(rstage, p_stages)
            p = spm12w_getp;
        end
        if ismember(rstage, [glm_stages,con_stages,rfx_stages])
            glm = spm12w_getp('type','glm');
        end
        if ismember(rstage, roi_stages)
            roi = spm12w_getp('type','roi');
        end
        if ismember(rstage, design_stages)
            des = spm12w_getp('type','des');
        end
    else
        % Loop through cell array of para files to load each one in
        for parfile = para_file
            spm12w_logger('msg',sprintf('[DEBUG] Loading parameters file: %s', ...
                                parfile{1}));
            % Figure out what type of file we're dealing with
            [~,paraname,~] = fileparts(parfile{1});
            ptype = strtok(paraname,'_');
            % A little bit of magic here since we don't know what parafile was given
            evalme = sprintf('%s = spm12w_getp(''type'',''%s'',''para_file'',''%s'');',...
                         ptype, ptype, parfile{1});
            eval(evalme);
        end
    end
    
    % Now run analysis stage in either parfor or regular loop
    % Note we cannot use subscripted references in parfor loops so we'll
    % use the usual indexing for both parfor and for loops (e.g., for i =
    % 1:length(x)
    
    % Enable to the parallel processing pool if necessary
    if gui.parloop == 1 && ismember(rstage, [p_stages,glm_stages,con_stages])
        guipool = parpool('local',gui.parcores);
    end
    
    % For preprocessing
    if ismember(rstage, p_stages)
        if gui.parloop == 1
            spm12w_logger('msg',sprintf(['[DEBUG] Starting parallel ',...
                          'processing on %d cores...'],gui.parcores))
            parfor sid_i = 1:length(sids)
                spm12w_preprocess('sid',sids{sid_i},'para_file',p.para_file);
            end
        else
            for sid_i = 1:length(sids)
                spm12w_preprocess('sid',sids{sid_i},'para_file',p.para_file);
            end
        end   
    end
    
    % For glm
    if ismember(rstage, glm_stages)
        if gui.parloop == 1
            spm12w_logger('msg',sprintf(['[DEBUG] Starting parallel ',...
                          'processing on %d cores...'],gui.parcores))
            parfor sid_i = 1:length(sids)
                spm12w_glm_compute('sid',sids{sid_i},'glm_file',glm.para_file);
            end
        else
            for sid_i = 1:length(sids)
                spm12w_glm_compute('sid',sids{sid_i},'glm_file',glm.para_file);
            end
        end   
    end
    
    % For contrasts
    if ismember(rstage, con_stages)
        if gui.parloop == 1
            spm12w_logger('msg',sprintf(['[DEBUG] Starting parallel ',...
                          'processing on %d cores...'],gui.parcores))
            parfor sid_i = 1:length(sids)
                spm12w_glm_contrast('sid',sids{sid_i},'glm_file',glm.para_file);
            end
        else
            for sid_i = 1:length(sids)
                spm12w_glm_contrast('sid',sids{sid_i},'glm_file',glm.para_file);
            end
        end   
    end
    
    % Close the parallel processing pool if it was opened
    if gui.parloop == 1 && ismember(rstage, [p_stages,glm_stages,con_stages])
        delete(guipool)
    end
    
    % For everything else there is no parfor or for looping over sids
    switch rstage     
        case 'rfx'
            spm12w_glm_rfx('sids',sids,'glm_file',glm.para_file);       

        case 'roi'
            spm12w_roitool('sids',sids,'roi_file',roi.para_file);    
        
        case 'dsearch'
            spm12w_designtool('type','search','des_file',des.para_file);
        
        case 'dcheck'
            spm12w_designtool('type','check','des_file',des.para_file);
        
        case 'dbuild'
            spm12w_designtool('type','build','des_file',des.para_file);
               
    end
    spm12w_logger('msg',gui.niceline)
    spm12w_logger('msg',sprintf('spm12w analysis completed on stage: %s', rstage))
return
