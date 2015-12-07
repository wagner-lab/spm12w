function spm12w(varargin)
% spm12w(stage, subjects)
%
% Inputs
% ------
% stage: String for the specific analysis stage desired. Supplying a stage
%        will override the gui and go directly to performing the analysis
%        stage.
%
% sids: Cell array of subject ids. If left unspecified, a dialog box will
%       appear asking the user to select subjects. If stage is empty, then
%       specifying sids will not override the gui. You may also specify
%       the path to a .mat file containing a saved sidlist cell array 
%       (e.g. sidlist = {'s01','s04','s12'}; save('sidlist.mat','sidlist'))
%
% Examples:
% Please note that spm12w should always be run from a study's top-level 
% directory.
%
%
% To select analyses via the spm12w gui, type spm12w with no argument:
%    
%       >> spm12w
%
% To run spm12w preprocessing on sids s01 and s02 overriding the gui:
%  
%       >> spm12w('stage','prep','sids',{'s01','s02'})
%
% To run spm12w preprocessing on sids s01 and s02 w/o overriding the gui:
%  
%       >> spm12w('sids',{'s01','s02'})
%
% To run spm12w preprocessing on a saved list of sids overriding the gui:
%  
%       >> spm12w('stage','prep','sids','./auxil/sidlist.mat')
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: May, 2015 | Updated: December, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('stage','', 'sids','');
args = spm12w_args('nargs',0, 'defaults', args_defaults, 'arguments', varargin);

% Check for user supplied sidlist mat file
if length(args.sids) == 1 && any(strfind(args.sids{1},'.mat'))
    args.sids = load(args.sids{1});
end

% If user submitted an analysis stage then they do not invoke gui.
if isempty(args.stage)
    ressurect_gui(args.sids) % Arise you shinning GUI.
else
    stager(args.stage, args.sids) 
end

function ressurect_gui(sids)
    close all
    % Define gui defaults
    % Define color pallet (push this to defaults later)
    gui.p_bcknd = [75,87,77] ./255;
    gui.p_panel = [141,152,118] ./255;
    gui.p_prepro = [175,66,10] ./255;
    
    % Define misc gui details
    gui.name = 'spm12w-gui';
    gui.imgfile = which('spm12w_title.jpg');
    % Define gui positions
    gui.pos_gui = [50,400,552,616]; %Position to place gui figure (X,Y,W,H)
    
    other = [203,179,69];
    other2 = [96,159,128];
    % Generate gui
    spm12wfig = figure('Name',gui.name, ...
                      'MenuBar','none', ...
                      'Toolbar','none', ...
                      'NumberTitle','off', ...
                      'Resize', 'off', ...
                      'Position',gui.pos_gui, ...
                      'color', gui.p_bcknd); 
    % Set gui image
    titleimage = axes('Units','pixels','position',[2,536,550,107]); 
    image(imread(gui.imgfile)); 
    axis('off','image') 
    % Define Panels
    ypos = 500; %Sets starting Y position. All other is relative to this
    % Main panel
    p0 = uipanel('BackgroundColor',(gui.p_panel),'Units','pixels',...
                 'Position',[52,1,450,568]); %p0: main panel
    % Preprocessing Panel
    p1 = uipanel('BackgroundColor',(gui.p_prepro),'Units','pixels',...
                 'Position',[60,ypos,435,62]); %p1: prepro panel
    % GLM / RFX Panel
    p2 = uipanel('BackgroundColor',([92,34,26] ./255),'Units','pixels',...
                 'Position',[60,(ypos-66),435,62]);
    % PPI Panel
    p3 = uipanel('BackgroundColor',([102,38,25] ./255),'Units','pixels',...
                 'Position',[60, (ypos-132),435,62]);
    % ROI and Conjunction panel
    p4 = uipanel('BackgroundColor',([112,42,24] ./255),'Units','pixels',...
                 'Position',[60, (ypos - 198),435,62]);
    % Segment and Dartel panel
    p5 = uipanel('BackgroundColor',([122,46,23] ./255),'Units','pixels',...
                 'Position',[60, (ypos - 236),435,34]);
    % REST preprocessing panel (wip)
    p6 = uipanel('BackgroundColor',([122,46,23] ./255),'Units','pixels',...
                 'Position',[60, (ypos - 302),435,62]);
    % Artifact detection panel (wip)
    p7 = uipanel('BackgroundColor',([112,42,24] ./255),'Units','pixels',...
                 'Position',[60, (ypos - 340),435,34]);   
    % Alphasim panel (wip)
    p8 = uipanel('BackgroundColor',([132,50,22] ./255),'Units','pixels',...
                 'Position',[60, (ypos - 378),435,34]);
    % Design search and build panel (wip)
    p9 = uipanel('BackgroundColor',([142,54,21] ./255),'Units','pixels',...
                 'Position',[60, (ypos - 416),435,34]);
    % Utilities panel
    p10 = uipanel('BackgroundColor',([152,58,20] ./255),'Units','pixels',...
                  'Position',[60, (ypos - 482),435,62]);
    % Define stage text 
    stage={'Preprocess','prep';
           'Preprocess & GLM','prep_glm';
           'Preprocess & GLM & Contrasts','prep_glm_con';
           'Estimate GLM','glm';
           'Estimate GLM & Contrasts','glm_con';
           'Compute Contrasts','con';
           'Random Effects','rfx';
           'VOI Extraction','voi';
           'WIP-PPI Reg Maker','ppi';
           'WIP-PPI Plotter','ppplot';
           'ROI Analysis','roi';
           'WIP-Conjunction Analysis','conj';
           'Segment Anatomy','seg8';
           'WIP-Template & Flow Fields','tempflow';
           'Normalize Anatomy','nanat';
           'Normalize Functional','nfunc';
           'WIP-REST Preprocess','rest';
           'ART Outlier Detection','art';
           'WIP-Monte Carlo Simulations (AlphaSim)','alphasim';
           'WIP-Design Search','dsearch';
           'WIP-Design Check','dcheck';
           'WIP-Design Build','dbuild';
           'WIP-Onset Maker','onsmk';
           'PAR/REC Converter', 'prepare'};
     %Set Buttons
     %Panel 01-Preprocessing
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{1,1},'Units', 'pixels','Position',[44,31,160,26],'Parent',p1,'Callback',{@stager, 'prep',sids});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{2,1},'Units', 'pixels','Position',[229,31,160,26],'Parent',p1,'Callback',{@stager, stage{2,2}});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{3,1},'Units', 'pixels','Position',[84,3,260,26],'Parent',p1,'Callback',{@stager, stage{3,2}});
     %Panel 02-GLM
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{4,1},'Units', 'pixels','Position',[44,31,160,26],'Parent',p2, 'Callback',{@stager, 4});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{5,1},'Units', 'pixels','Position',[229,31,160,26],'Parent',p2, 'Callback',{@stager, 5});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{6,1},'Units', 'pixels','Position',[44,3,160,26],'Parent',p2, 'Callback',{@choice, 6});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{7,1},'Units', 'pixels','Position',[229,3,160,26],'Parent',p2,'Callback',{@choice, 7});
     %Panel 03-PPI
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{8,1},'Units', 'pixels','Position',[44,31,160,26],'Parent',p3, 'Callback',{@choice, 8});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{9,1},'Units', 'pixels','Position',[229,31,160,26],'Parent',p3, 'Callback',{@choice, 9});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{10,1},'Units', 'pixels','Position',[84,3,260,26],'Parent',p3, 'Callback',{@choice, 10});      
     %Panel 04-ROI and Conjunction
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{11,1},'Units', 'pixels','Position',[44,3,160,26],'Parent',p4, 'Callback',{@choice, 11});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{12,1},'Units', 'pixels','Position',[229,3,160,26],'Parent',p4,'Callback',{@choice, 12});
     %Panel 05-Segment and DARTEL
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{13,1},'Units', 'pixels','Position',[44,31,160,26],'Parent',p5, 'Callback',{@choice, 13});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{14,1},'Units', 'pixels','Position',[229,31,160,26],'Parent',p5, 'Callback',{@choice, 14});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{15,1},'Units', 'pixels','Position',[44,3,160,26],'Parent',p5, 'Callback',{@choice, 15});      
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{16,1},'Units', 'pixels','Position',[229,3,160,26],'Parent',p5, 'Callback',{@choice, 16});           
     %Panel 06-REST Prestageing
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{17,1},'Units', 'pixels','Position',[84,3,260,26],'Parent',p6,'Callback',{@choice, 17});
     %Panel 07-Artifact Detection
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{18,1},'Units', 'pixels','Position',[84,3,260,26],'Parent',p7,'Callback',{@choice, 18});
     %Panel 08-AlphaSim
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{19,1},'Units', 'pixels','Position',[84,3,260,26],'Parent',p8,'Callback',{@choice, 19});
     %Panel 09-Design Search, Check & Build
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{20,1},'Units', 'pixels','Position',[44,31,160,26],'Parent',p9,'Callback',{@choice, 20});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{21,1},'Units', 'pixels','Position',[229,31,160,26],'Parent',p9,'Callback',{@choice, 21});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{22,1},'Units', 'pixels','Position',[84,3,260,26],'Parent',p9,'Callback',{@choice, 22});
     %Panel 10-Utilities
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{23,1},'Units', 'pixels','Position',[44,3,160,26],'Parent',p10, 'Callback',{@choice, 23});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',stage{24,1},'Units', 'pixels','Position',[229,3,160,26],'Parent',p10, 'Callback',{@choice, 24});
return

% Implement the requested stage
function stager(h, eventdata, rstage, sids)
    close(gcf);

    % Load sids for stages that require them
    sidstages = {'prep','prep_glm','prep_glm_con','glm','glm_con','con',...
                 'rfx','voi','roi','seg8','nanat','nfunc','rest'};
    if isempty(sids)
        if ismember(rstage,sidstages(1:3))
            sids = spm12w_getsid(fullfile(pwd,'raw'));
        elseif ismember(rstage,sidstages(4:end))
            sids = spm12w_getsid();
        end
    end
    
    switch rstage
        case 'prep'
            p = spm12w_getp;
            for sid = sids
                spm12w_preprocess('sid',sid{1},'para_file',p.para_file);
            end

        case 'prep_glm'
            p = spm12w_getp('type','p');
            g = spm12w_getp('type','glm');
            for sid = sids
                spm12w_preprocess('sid',sid{1},'para_file',p.para_file);
                spm12w_glm_compute('sid',sid{1},'glm_file',g.para_file);
            end

        case 'prep_glm_con'
            p = spm12w_getp;
            g = spm12w_getp('type','glm');
            for sid = sids
                spm12w_preprocess('sid',sid{1},'para_file',p.para_file);
                spm12w_glm_compute('sid',sid{1},'glm_file',g.para_file);
                spm12w_contrasts('sid',sid{1},'glm_file',g.para_file);
            end

        case 'glm'
            g = spm12w_getp('type','glm');
            for sid = sids
                spm12w_glm_compute('sid',sid{1},'glm_file',g.para_file);
            end

        case 'glm_con'
            g = spm12w_getp('type','glm');
            for sid = sids
                spm12w_glm_compute('sid',sid{1},'glm_file',g.para_file);
                spm12w_contrasts('sid',sid{1},'glm_file',g.para_file);
            end
        case 'con'
            g = spm12w_getp('type','glm');
            for sid = sids
                spm12w_contrasts('sid',sid{1},'glm_file',g.para_file);
            end

        case 'rfx'
            p = spm12w_getp;
            for sid = sids
                spm12w_preprocess('sid',sid{1},'para_file',p.para_file);
            end

        case 'roi'
            p = spm12w_getp;
            for sid = sids
                spm12w_preprocess('sid',sid{1},'para_file',p.para_file);
            end
        
        otherwise
            error('not supported')        
    end
return
