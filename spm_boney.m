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

global boney

% Test for CAT installation 
catdir = fullfile( spm('dir') , 'toolbox' , 'cat12'); 
if ~exist( catdir , 'dir' )
  error('Error:NoCAT:(',sprintf('Cannot see CAT12 directory: \n%s', catdir))
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
  if ~isfield(boney,'expertgui') 
    boney.expertgui = 0; 
  end
end

if restart
  spm 'fmri'; spm_boney;
end
  
addpath(fullfile(spm('dir'),'toolbox','boney'));
rev = '0.1';

mode = {' ',' Expert mode',' Developer mode'}; 
SPMid = spm('FnBanner',mfilename,rev);
spm('FnUIsetup',['BONEY' mode{boney.expertgui + 1}]); 
F   = spm_figure('GetWin'); spm_figure('clear',F); 
ax  = axes(F,'Position',[0 0 1 1]); 
img = imread( fullfile( spm('Dir'), 'toolbox', 'boney', 'docs', 'Kalc-HBM2023-Scull.jpg') ); 
imshow(img); 


% main menu
% -------------------------------------------------------------------------
fig = spm_figure('GetWin','Interactive');
h0  = uimenu(fig,...
  'Label',            'Boney',...
  'Separator',        'on',...
  'Tag',              'boney',...
  'HandleVisibility', 'on');

h1 = uimenu(h0,...
  'Label',            'Bone segmentation',...
  'Separator',        'off',...
  'Tag',              'Bone segmentation',...
  'CallBack',         'spm_jobman(''interactive'','''',''spm.tools.boney.segment'');',...
  'HandleVisibility', 'on');

h2 = uimenu(h0,...
  'Label',            'Export results as CSV',...
  'Separator',        'off',...
  'Tag',              'CSV export',...
  'CallBack',         'spm_jobman(''interactive'','''',''spm.tools.cat.tools.xml2csv'','''');', ...
  'HandleVisibility', 'off');

h3 = uimenu(h0,...
  'Label',            'Show help',...
  'Separator',        'on',...
  'Tag',              'Help',...
  'CallBack',         ['' ...
    'F   = spm_figure(''GetWin''); ' ...
    'spm_figure(''clear'',F); ' ...
    'JW  = findobj(F.Children,''Type'',''hgjavacomponent''); ' ...
    'JW.Visible = ''off''; ' ...
    '' ...
    'ax  = axes(F,''Position'',[0 0 1 1]); ' ...
    'img = imread( fullfile( spm(''Dir''), ''toolbox'', ''boney'', ''docs'', ''Kalc-HBM2023-Scull.jpg'') ); ' ...
    'imshow(img); '], ...
  'HandleVisibility', 'off');
 

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
end
