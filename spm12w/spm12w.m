function spm12w(varargin)
% spm12w(stage, subjects)
%
% Inputs
% ------
% stage: Number (e.g., 3) or vector of numbers (e.g. [1,2,3,4,5]) for the 
%        specific analysis stage desired. Supplying a stage will override
%        the gui and go directly to performing the analysis stage.
%
% sids: Cell array of subject ids. If left unspecified, a dialog box will
%       appear asking the user to select subjects. If stage is empty, then
%       specifying sids will not override the gui. 
%

% Examples:
% Please note that spm12w should always be run from a study's top-level 
% directory.
%
% To run spm12w preprocessing on subjects s01 and s02 overriding the gui:
%  
%       >> spm12w('stage',1,'sids',{'s01','s02'})
%
% To select analyses via the spm12w gui, type spm12w with no argument:
%    
%       >> spm12w
%
% # spm12w was developed by the Wagner, Heatherton & Kelley Labs
% # Author: Dylan Wagner | Created: May, 2015 | Updated: December, 2015
% =======1=========2=========3=========4=========5=========6=========7=========8

% Parse inputs
args_defaults = struct('stage','', 'sids','');
args = spm12w_args('nargs',0, 'defaults', args_defaults, 'arguments', varargin);

% If user submitted an analysis stage then they do not invoke gui.
if isempty(args.stage)
    spm12w_stager(args.stage, args.sids)
else
    % Arise you shining GUI! Arise!
    ressurect_gui()
end

function ressurect_gui()
    close all
    % Generate figure
    spm12wfig = figure('Name','spm12w-gui', ...
                      'MenuBar','none', ...
                      'Toolbar','none', ...
                      'NumberTitle','off', ...
                      'Resize', 'off', ...
                      'Position',[50,400,552,616]); %X,Y,W,H
    %Change background
    set(spm12wfig, 'Color', ([202,198,191] ./255));
    %Set image
    imgfile    = which('spm12w_title.jpg');
    titleimage = axes('Units', 'pixels','position', [2,536,550,107]); %X,Y then Width, Height
    titleimg   = imread(imgfile);
    image(titleimg); 
    axis off 
    axis image
    %Set Panels
    ypos = 500; %Sets starting Y position of panels. All other position is relative
    p0 = uipanel('BackgroundColor',([47,48,48] ./255),'Units','pixels','Position',[52, 1,450,568]);
    p1 = uipanel('BackgroundColor',([82,30,27] ./255),'Units','pixels','Position',[60, ypos,435,62]);
    p2 = uipanel('BackgroundColor',([92,34,26] ./255),'Units','pixels','Position',[60, (ypos - 66),435,62]);
    p3 = uipanel('BackgroundColor',([102,38,25] ./255),'Units','pixels','Position',[60, (ypos - 104),435,34]);
    p4 = uipanel('BackgroundColor',([112,42,24] ./255),'Units','pixels','Position',[60, (ypos - 170),435,62]);
    p5 = uipanel('BackgroundColor',([122,46,23] ./255),'Units','pixels','Position',[60, (ypos - 208),435,34]);
    p6 = uipanel('BackgroundColor',([122,46,23] ./255),'Units','pixels','Position',[60, (ypos - 274),435,62]);
    p7 = uipanel('BackgroundColor',([112,42,24] ./255),'Units','pixels','Position',[60, (ypos - 312),435,34]);   
    p8 = uipanel('BackgroundColor',([132,50,22] ./255),'Units','pixels','Position',[60, (ypos - 350),435,34]);
    p9 = uipanel('BackgroundColor',([142,54,21] ./255),'Units','pixels','Position',[60, (ypos - 388),435,34]);
    p10 = uipanel('BackgroundColor',([152,58,20] ./255),'Units','pixels','Position',[60, (ypos - 454),435,62]);
    p11 = uipanel('BackgroundColor',([162,64,19] ./255),'Units','pixels','Position',[60, (ypos - 492),435,34]);   
    
    process={'Preprocess',...                         %1
             'Preprocess/GLM',...                     %2
             'Preprocess/GLM/Contrasts',...           %3
             'Estimate GLM',...                       %4
             'Estimate GLM/Contrasts',...             %5
             'Compute Contrasts',...                  %6
             'Random Effects',...                     %7
             'VOI Extraction',...                     %8
             'PPI Reg Maker',...                      %9
             'PPI Plotter',...                        %10
             'ROI Analysis',...                       %11
             'Conjunction Analysis',...               %12
             'Segment Anatomy',...                    %13
             'Template & Flow Fields',...             %14
             'Normalize Anatomy',...                  %15
             'Normalize Functional',...               %16
             'REST Preprocess',...                    %17
             'ART Outlier Detection',...              %18
             'Monte Carlo Simulations (AlphaSim)',... %19
             'Design Search',...                      %20
             'Design Check',...                       %21
             'Design Build',...                       %22
             'Onset Maker',...                        %23
             'PAR/REC Converter'};                    %24
     %Set Buttons
     %Panel 01-Preprocessing
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{1},'Units', 'pixels','Position',[44,31,160,26],'Parent',p1,'Callback',{@choice, 1});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{2},'Units', 'pixels','Position',[229,31,160,26],'Parent',p1,'Callback',{@choice, 2});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{3},'Units', 'pixels','Position',[84,3,260,26],'Parent',p1,'Callback',{@choice, 3});
     %Panel 02-GLM
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{4},'Units', 'pixels','Position',[44,31,160,26],'Parent',p2, 'Callback',{@choice, 4});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{5},'Units', 'pixels','Position',[229,31,160,26],'Parent',p2, 'Callback',{@choice, 5});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{6},'Units', 'pixels','Position',[84,3,260,26],'Parent',p2, 'Callback',{@choice, 6});
     %Panel 03-RFX
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{7},'Units', 'pixels','Position',[84,3,260,26],'Parent',p3,'Callback',{@choice, 7});
     %Panel 04-PPI
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{8},'Units', 'pixels','Position',[44,31,160,26],'Parent',p4, 'Callback',{@choice, 8});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{9},'Units', 'pixels','Position',[229,31,160,26],'Parent',p4, 'Callback',{@choice, 9});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{10},'Units', 'pixels','Position',[84,3,260,26],'Parent',p4, 'Callback',{@choice, 10});      
     %Panel 05-ROI and Conjunction
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{11},'Units', 'pixels','Position',[44,3,160,26],'Parent',p5, 'Callback',{@choice, 11});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{12},'Units', 'pixels','Position',[229,3,160,26],'Parent',p5,'Callback',{@choice, 12});
     %Panel 06-Segment and DARTEL
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{13},'Units', 'pixels','Position',[44,31,160,26],'Parent',p6, 'Callback',{@choice, 13});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{14},'Units', 'pixels','Position',[229,31,160,26],'Parent',p6, 'Callback',{@choice, 14});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{15},'Units', 'pixels','Position',[44,3,160,26],'Parent',p6, 'Callback',{@choice, 15});      
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{16},'Units', 'pixels','Position',[229,3,160,26],'Parent',p6, 'Callback',{@choice, 16});           
     %Panel 07-REST Preprocessing
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{17},'Units', 'pixels','Position',[84,3,260,26],'Parent',p7,'Callback',{@choice, 17});
     %Panel 08-Artifact Detection
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{18},'Units', 'pixels','Position',[84,3,260,26],'Parent',p8,'Callback',{@choice, 18});
     %Panel 09-AlphaSim
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{19},'Units', 'pixels','Position',[84,3,260,26],'Parent',p9,'Callback',{@choice, 19});
     %Panel 10-Design Search, Check & Build
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{20},'Units', 'pixels','Position',[44,31,160,26],'Parent',p10,'Callback',{@choice, 20});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{21},'Units', 'pixels','Position',[229,31,160,26],'Parent',p10,'Callback',{@choice, 21});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{22},'Units', 'pixels','Position',[84,3,260,26],'Parent',p10,'Callback',{@choice, 22});
     %Panel 11-Utilities
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{23},'Units', 'pixels','Position',[44,3,160,26],'Parent',p11, 'Callback',{@choice, 23});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{24},'Units', 'pixels','Position',[229,3,160,26],'Parent',p11, 'Callback',{@choice, 24});
return


