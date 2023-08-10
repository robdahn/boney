function spm_boney(expertgui)
% Toolbox wrapper to call boney functions
% _________________________________________________________________________
%
% Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________

% TODO: 
% * feature:  add soft skull-crossed-bones image in the progress window 
% * feature:  add default parameter file? eg. to highligh changes from default?
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
  
  Pposter = fullfile( spm('Dir'), 'toolbox', 'boney', 'docs', 'Kalc-HBM2023-Scull.jpg'); 
  Phelp   = fullfile( spm('Dir'), 'toolbox', 'boney', 'README.html'); 
  mode = {' ',' Expert mode',' Developer mode'}; 
  SPMid = spm('FnBanner',mfilename,rev);
  spm('FnUIsetup',['BONEY' mode{boned.expertgui + 1}]); 
  % just add some image for the SPM figure that is not suporting the html 
  % help in new Matlab version - assure the figure size!
  F = spm_figure('GetWin'); 
  spm_figure('clear',F); 
  Fpos = get(F,'Position'); 
  h = imshow(imread(Pposter)); 
  set(get(h,'Parent'),'Position',[0 0 1 1]);
  set(F,'Position',Fpos);
  
  % main menu
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
      ['The GUI may not working under Octave yet. To create a result table call: \n\n '...
       '   spm_jobman(''interactive'','''',''spm.tools.boney.xml2csv''); \n\n']); 
  end

end
