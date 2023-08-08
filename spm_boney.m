function spm_boney(expertgui)
% Toolbox wrapper to call boney functions
% _________________________________________________________________________
%
% Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________

%#ok<*NASGU>

global boney %#ok<GVMIS> 

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
  if isfield(boney,'expertgui') && boney.expertgui ~= expertgui
    restart = 1; 
  else
    restart = 0;
  end
  boney.expertgui = expertgui; 
else
  restart = 0;
  boney.expertgui = 0; 
end

if restart
  spm 'fmri'; spm_boney;
end
  
addpath(fullfile(spm('dir'),'toolbox','boney'));
rev = '0.1';

Pposter = fullfile( spm('Dir'), 'toolbox', 'boney', 'docs', 'Kalc-HBM2023-Scull.jpg'); 
Phelp   = fullfile( spm('Dir'), 'toolbox', 'boney', 'README.html'); 
mode = {' ',' Expert mode',' Developer mode'}; 
SPMid = spm('FnBanner',mfilename,rev);
spm('FnUIsetup',['BONEY' mode{boney.expertgui + 1}]); 
F = spm_figure('GetWin'); spm_figure('clear',F); 
h = imshow(Pposter); 
set(get(h,'Parent'),'Position',[0 0 1 1]);

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
  'CallBack',         'spm_jobman(''interactive'','''',''spm.tools.cat.tools.xml2csv'','''');', ...
  'HandleVisibility', 'off');

h3 = uimenu(h0,...
  'Label',            'Display intro',...
  'Separator',        'on',...
  'Tag',              'Intro',...
  'CallBack',         ['' ...
    'F   = spm_figure(''GetWin''); ' ...
    'spm_figure(''clear'',F); ' ...
    'JW  = findobj(F.Children,''Type'',''hgjavacomponent''); ' ...
    'JW.Visible = ''off''; ' ...
    'imshow(''' Pposter '''); ' ...
    'set(get(h,''Parent''),''Position'',[0 0 1 1]); ' ], ...
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
if boney.expertgui==1
  cat_io_cprintf([0.0 0.0 1.0],sprintf('(Expert mode)')); 
elseif boney.expertgui==2
  cat_io_cprintf([1.0 0.0 0.0],sprintf('(Developer mode)')); 
end
fprintf('\n\n');
 
if strcmpi(spm_check_version,'octave') 
% in Ocatve the menu is not working  
  spm_jobman('interactive','','spm.tools.boney.segment');
  fprintf('warn',...
    ['The GUI may not working under Octave yet. To create a result table call: \n\n '...
     '   spm_jobman(''interactive'','''',''spm.tools.cat.tools.xml2csv''); \n\n']); 
end

end
