function spm_boney(expertgui)
% Toolbox wrapper to call boney functions
% _________________________________________________________________________
% _________________________________________________________________________
% Robert Dahnke, 20191003

%#ok<*NASGU>

global boney

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
[~,Fgraph,~] = spm('FnUIsetup',['BONEY' mode{boney.expertgui + 1}]);
url = fullfile(spm('Dir'),'toolbox','boney','html','boney.html');
spm_help('!Disp',url,'',Fgraph,'Boney toolbox for SPM');



% main menu
fig = spm_figure('GetWin','Interactive');
h0  = uimenu(fig,...
  'Label',            'boney',...
  'Separator',        'on',...
  'Tag',              'boney',...
  'HandleVisibility', 'on');


% volume
% -------------------------------------------------------------------------
h1 = uimenu(h0,...
  'Label',            'Bone segmentation',...
  'Separator',        'off',...
  'Tag',              'Bone segmentation',...
  'CallBack',         'spm_jobman(''interactive'','''',''spm.tools.boney.segment'');',...
  'HandleVisibility', 'on');



% manual/help
% -------------------------------------------------------------------------
h2  = uimenu(h0,...
  'Label',            'Open help',...
  'Separator',        'off',...
  'Tag',              'Help',...
  'CallBack',         'spm_help(''!Disp'',fullfile(spm(''Dir''),''toolbox'',''boney'',''boney.html''),'''',spm_figure(''GetWin''),''Boney toolbox for SPM'');', ...
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

