function boney = tbx_cfg_boney(expertgui)
% Configuration file for boney toolbox batches
% _________________________________________________________________________
% _________________________________________________________________________
% Robert Dahnke 2023
%
%#ok<*AGROW,*INUSD>

if nargin==0, expertgui = 1; end

addpath(fileparts(which(mfilename)));

% try to estimate number of processor cores
try
  numcores = cat_get_defaults('extopts.nproc');
  % because of poor memory management use only half of the cores for windows
  if ispc
    numcores = round(numcores/2);
  end
  numcores = max(numcores,1);
catch
  numcores = 0;
end

% force running in the foreground if only one processor was found or for compiled version
% or for Octave
if numcores == 1 || isdeployed || strcmpi(spm_check_version,'octave'), numcores = 0; end

% parallelize
% ____________________________________________________________________
nproc         = cfg_entry;
nproc.tag     = 'nproc';
nproc.name    = 'Split job into separate processes';
nproc.strtype = 'w';
nproc.val     = {numcores};
nproc.num     = [1 1];
nproc.help    = {
    'In order to use multi-threading the CAT12 segmentation job with multiple subjects can be split into separate processes that run in the background. You can even close Matlab, which will not affect the processes that will run in the background without GUI. If you do not want to run processes in the background then set this value to 0.'
    ''
    'Please note that no additional modules in the batch can be run except CAT12 segmentation. Any dependencies will be broken for subsequent modules.'
  };


%  basic entries
%  ------------------------------------------------------------------------

files              = cfg_files;
files.tag          = 'files';
files.name         = 'Images';
files.help         = {'Select images that should be processed.'};
files.filter       = 'image';
files.ufilter      = '.*';
files.num          = [1 Inf];

% verb 
verb                = cfg_menu;
verb.tag            = 'verb';
verb.name           = 'Verbose processing';
verb.labels         = {'Yes','No'}; 
verb.values         = {1,0};
verb.val            = {1}; 
verb.help           = {
  'Verbose processing. '
  ''};


%  batches
%  ------------------------------------------------------------------------
segment = boney_cfg_segment(files,nproc,expertgui,verb);

%  main
%  ------------------------------------------------------------------------
boney               = cfg_choice;
boney.name          = 'Boney';
boney.tag           = 'boney';
boney.values        = {segment};



%  subfunctions
%  ------------------------------------------------------------------------

function segment = boney_cfg_segment(files,nproc,expertgui,verb)

% opts.
% .method 
% .res
% .msk
% .atlas
% 
% output
% .images
% .surfaces


% == opts == 
  % development option - start with SPM, will add CAT that is maybe more
  % robust in some cases, maybe use different presettings, eg TPMs for
  % children later 
  pmethod               = cfg_menu;
  pmethod.tag           = 'pmethod';
  pmethod.name          = 'Preprocessing method';
  pmethod.labels        = {'SPM','CAT'};
  pmethod.values        = {'spm','cat'};
  pmethod.val           = {'spm'};
  pmethod.hidden        = expertgui<2;
  pmethod.help          = {'Preprocessing method to segment the different tissue classes in the given image.';''};

  % main method of this toolbox
  bmethod               = cfg_menu;
  bmethod.tag           = 'bmethod';
  bmethod.name          = 'Bone precessing method';
  bmethod.labels        = {'Volume-based','Surface-based'};
  bmethod.values        = {'med','cat'};
  if expertgui
    bmethod.labels      = [ { 'SPM (first prototype)' , 'Volume-based (first prototype)' } , bmethod.labels ];
    bmethod.values      = [ { 'SPM0' , 'MED0' ,  'cat'} , bmethod.values ];
  end
  bmethod.val           = {'spm'};
  bmethod.hidden        = expertgui<1;
  bmethod.help          = {'Bone processing method.';''};

  affreg               = cfg_menu;
  affreg.tag           = 'affreg';
  affreg.name          = 'Apply affine registration';
  affreg.labels        = {'No','Yes'};
  affreg.values        = {0,1};
  affreg.val           = {0};
  affreg.hidden        = expertgui<2;
  affreg.help          = {'Apply bone-focused affine registration.';''};

  reduce               = cfg_menu;
  reduce.tag           = 'reduce';
  reduce.name          = 'Processing resolution';
  reduce.labels        = {'full','half','quater'};
  reduce.values        = {1,2,4};
  reduce.val           = {2};
  reduce.hidden        = expertgui<1;
  reduce.help          = {'Processing resolution.';''};

  mask               = cfg_menu;
  mask.tag           = 'Pmask';
  mask.name          = 'Bone mask';
  mask.help          = {'Bone atlas that describe the major skull bones (frontal, pariatal-left/right, temporal, occipital). '};
  mask.filter        = {'image'};
  mask.dir           = fullfile(spm('dir'),'toolbox','boney','boney_KADA_bone-mask.nii'); 
  mask.ufilter       = '.*';
  mask.num           = [0 1];
  mask.hidden        = expertgui<1;
  mask.help          = {'Mask problematic regions.';''};
 
  atlas              = cfg_files;
  atlas.tag          = 'Patlas';
  atlas.name         = 'Bone atlas';
  atlas.help         = {'Bone atlas that describe the major skull bones (frontal, pariatal-left/right, temporal, occipital). '};
  atlas.filter       = {'image'};
  atlas.dir          = fullfile(spm('dir'),'toolbox','boney','boney_KADA_bone-regions.nii'); 
  atlas.ufilter      = '.*';
  atlas.num          = [0 1];
  atlas.hidden       = expertgui<1;
  %atlas.val          = {{''}};

  % main parameter structure passed to cat_plot_boxplot
  opts           = cfg_exbranch;
  opts.tag       = 'opts';
  opts.name      = 'Options';
  opts.val       = { pmethod , bmethod, nproc; affreg; verb }; 
  opts.help      = {'Specify processing parameters.' ''};



% == output ==  
  report               = cfg_menu;
  report.tag           = 'report';
  report.name          = 'Write report';
  report.labels        = {'No','Values','Volumes','Surfaces'};
  report.values        = {0,1,2,3};
  report.val           = {3};
  report.hidden        = expertgui<1;
  report.help          = {'Write JPG report file with values only (fast), additional volume or surface overlays if available.';''};

  writevol               = cfg_menu;
  writevol.tag           = 'writevol';
  writevol.name          = 'Write volumes';
  writevol.labels        = {'No','Yes'};
  writevol.values        = {0,1};
  writevol.val           = {0};
  writevol.hidden        = expertgui<1;
  writevol.help          = {'Write refined bone segment volumes used for extraction.';''};
 
  writesurf               = cfg_menu;
  writesurf.tag           = 'writesurf';
  writesurf.name          = 'Write surfaces';
  writesurf.labels        = {'No','Yes'};
  writesurf.values        = {0,1};
  writesurf.val           = {0};
  writesurf.hidden        = expertgui<1;
  writesurf.help          = {'Write bone surfaces used for extraction (if processed).';''};

  output           = cfg_exbranch;
  output.tag       = 'output';
  output.name      = 'Output';
  output.val       = { report , writevol , writesurf }; 
  output.help      = {'Specify output parameters.' ''};




  segment          = cfg_exbranch;
  segment.tag      = 'segment';
  segment.name     = 'Bone segmentation';
  segment.val      = {files,opts,output}; 
  segment.prog     = @boney_segment;
  %segment.vout     = @vout_boney;
  segment.help     = {
   '' ''};


  