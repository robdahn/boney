function spm_boney(expertgui)
%spm_boney. Toolbox wrapper to call boney functions.
%
%  spm_boney(expertgui)
% 
%
% Creation of a global variable "boned" for default settings. 
% Setting of GUI level (0-default, 1-expert, 2-developer)
% _________________________________________________________________________
%
% Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________

% TODO: 
% * feature:  add default parameter file? e.g. to highligh changes from default?
%   - not realy required - keep it simple

%#ok<*NASGU>
  
  global boned %#ok<GVMIS> 
  
  % Test for CAT installation 
  catdir = fullfile( spm('dir') , 'toolbox' , 'cat12'); 
  if ~exist( catdir , 'dir' )
    error('Error:noCATtoolbox','Cannot see CAT12 directory: \n%s', catdir)
  end
  
  %% define default parameters
  if exist('expertgui','var')
    switch expertgui
      case {0,'default'},   expertgui = 0; 
      case {1,'expert'},    expertgui = 1; 
      case {2,'developer'}, expertgui = 2; 
      otherwise 
        expertgui = 0; 
    end
    if isfield(boned,'expertgui') && boned.expertgui ~= expertgui
      restart = 1; 
    else
      restart = 0;
    end
    boned.expertgui = expertgui; 
  else
    cat_io_cprintf([1 0 0],'\n  Currently we are all boney developer if not stated otherwise!\n\n');
    restart = 0;
    boned.expertgui = 2; 
  end
  
  if restart
    spm 'fmri'; spm_boney(boned.expertgui);
    return
  end
    
  addpath(fullfile(spm('dir'),'toolbox','boney'));
  rev = '0.1';
  
  % skull: Human skull drawing, medical vintage illustration psd. Free public domain CC0 image.
  Pposter = fullfile( spm('Dir'), 'toolbox', 'boney', 'docs', 'Kalc-HBM2023-Scull.jpg'); 
  Pinter  = fullfile( spm('Dir'), 'toolbox', 'boney', 'images', 'SPM_progress_skullsm.png'); 
  Phelp   = fullfile( spm('Dir'), 'toolbox', 'boney', 'README.html'); 
  mode = {' ',' Expert mode',' Developer mode'}; 
  SPMid = spm('FnBanner',mfilename,rev);
  spm('FnUIsetup',['BONEY' mode{boned.expertgui + 1}]); 
  % just add some image for the SPM figure that is not supporting the html 
  % help in new Matlab version - assure the figure size!
  F = spm_figure('GetWin'); 
  spm_figure('clear',F); 
  Fpos = get(F,'Position'); 
  h = imshow(imread(Pposter)); 
  set(get(h,'Parent'),'Position',[0 0 1 1]);
  set(F,'Position',Fpos);


  %% Progress Window figure
  Finter = spm_figure('GetWin','Interactive');
  spm_figure('clear',Finter); 
  [img, ~, alphachannel] = imread(Pinter);
  h = image(img, 'AlphaData', alphachannel); axis equal off
  set(get(h,'Parent'),'Position',[0 0.05 1 0.7]);
  % SPM watermark
  spm_figure('WaterMark',Finter,'Boney','WaterMark',1);
  WM = findobj('Tag','WaterMark'); 
  WM2 = copyobj(WM,get(WM,'Parent'));
  set(WM,'Position',[0.45 0.795 .1 .1 ]);
  set(WM2,'Position',[0.45 0.80 .1 .1 ]);
  WMF = get(WM,'Children'); 
  set(WMF,'Color',ones(1,3)*.5,'FontSize', 80.5); %,'FontName','Arial'); 
  WMF2 = get(WM2,'Children'); 
  set(WMF2,'Color',ones(1,3)*.95,'FontSize', 80); %,'FontName','Arial'); 


  %% main menu
  % -------------------------------------------------------------------------
  fig = spm_figure('GetWin','Interactive');
  h0  = uimenu(fig,...
    'Label',            'Boney',...
    'Separator',        'on',...
    'Tag',              'boney',...
    'HandleVisibility', 'on');
  
  h1 = uimenu(h0,...
    'Label',            'Bone/head preprocessing',...
    'Separator',        'off',...
    'Tag',              'Bone/head preprocessing',...
    'CallBack',         'spm_jobman(''interactive'','''',''spm.tools.boney.segment'');',...
    'HandleVisibility', 'on');
  
  h2 = uimenu(h0,...
    'Label',            'Export results as CSV',...
    'Separator',        'off',...
    'Tag',              'CSV export',...
    'CallBack',         'spm_jobman(''interactive'','''',''spm.tools.boney.xml2csv'','''');', ...
    'HandleVisibility', 'off');
  
  h3 = uimenu(h0,...
    'Label',            'Display intro',...
    'Separator',        'on',...
    'Tag',              'Intro',...
    'CallBack',         ['' ...
      'F   = spm_figure(''GetWin''); ' ...
      'spm_figure(''clear'',F); ' ...
      'Fpos = get(F,''Position''); ' ...
      'h = imshow(''' Pposter '''); ' ...
      'set(get(h,''Parent''),''Position'',[0 0 1 1]); ' ...
      'set(F,''Position'',Fpos);' ], ...
    'HandleVisibility', 'off');
  
  h4 = uimenu(h0,...
    'Label',            'Show help',...
    'Separator',        'off',...
    'Tag',              'Help',...
    'CallBack',         ['web(''' Phelp ''')']); 
   
  
  %% command line output
  % -------------------------------------------------------------------------
  cat_io_cprintf('silentreset');
  cat_io_cprintf([0.0 0.0 0.5],sprintf([ ...
      '    \n' ...
      '    Boney Toolbox\n' ...
      '    Version ' rev ' alpha ']));
  if boned.expertgui==1
    cat_io_cprintf([0.0 0.0 1.0],sprintf('(Expert mode)')); 
  elseif boned.expertgui==2
    cat_io_cprintf([1.0 0.0 0.0],sprintf('(Developer mode)')); 
  end
  fprintf('\n\n');
   
  if strcmpi(spm_check_version,'octave') 
  % in Ocatve the menu is not working  
    spm_jobman('interactive','','spm.tools.boney.segment');
    fprintf(...
      ['The GUI may not be working under Octave yet. To create a result table call: \n\n '...
       '   spm_jobman(''interactive'','''',''spm.tools.boney.xml2csv''); \n\n']); 
  end

end
