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
  
  % force running in the foreground if only one processor was found or for compiled version or for Octave
  if numcores == 1 || isdeployed || strcmpi(spm_check_version,'octave'), numcores = 0; end
  
  % parallelize
  % ____________________________________________________________________
  nproc         = cfg_entry;
  nproc.tag     = 'nproc';
  nproc.name    = 'Split job into separate processes';
  nproc.strtype = 'w';
  nproc.val     = {numcores};
  nproc.num     = [1 1];
  nproc.hidden  = true; 
  nproc.help    = {
     ['NOT YET IMPLEMENTED: ' ...
      'In order to use multi-threading jobs with multiple subjects can be split into separate processes that run in the background. ' ...
      'You can even close Matlab, which will not affect the processes that will run in the background without GUI. ' ...
      'If you do not want to run processes in the background then set this value to 0. ' ]
      ''
      'Please note that no additional modules in the batch can be run. Any dependencies will be broken for subsequent modules.'
    };
  
  
  %  basic entries
  %  ------------------------------------------------------------------------
  files               = cfg_files;
  files.tag           = 'files';
  files.name          = 'Images';
  files.help          = {'Select images that should be processed.'};
  files.filter        = 'image';
  files.ufilter       = '.*';
  files.num           = [1 Inf];
   
  verb                = cfg_menu;
  verb.tag            = 'verb';
  verb.name           = 'Verbose processing';
  verb.labels         = {'No','Yes'};
  verb.values         = {0,1};
  if expertgui 
    verb.labels       = [ verb.labels, {'Yes - Details'} ];
    verb.values       = [ verb.values, {2} ];
  end
  verb.val            = {1}; 
  verb.hidden         = expertgui<1; 
  
  
  %  batches
  %  ------------------------------------------------------------------------
  segment = boney_cfg_segment(files,nproc,expertgui,verb);
  
  
  %  main
  %  ------------------------------------------------------------------------
  boney               = cfg_choice;
  boney.name          = 'Boney';
  boney.tag           = 'boney';
  boney.values        = {segment};
end


%  subfunctions
%  ------------------------------------------------------------------------
function segment = boney_cfg_segment(files,nproc,expertgui,verb)
%boney_cfg_segment. Bone segmentation and value extraction.

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
  pmethod.hidden        = expertgui<2; %%%%%%
  pmethod.help          = {'Preprocessing method to segment the different tissue classes in the given image.';''};

  % main method of this toolbox
  bmethod               = cfg_menu;
  bmethod.tag           = 'bmethod';
  bmethod.name          = 'Bone processing method';
  bmethod.labels        = {'SPM mat-file', 'Volume-based', 'Surface-based'};
  bmethod.values        = {1, 2, 3};
  if expertgui
    bmethod.labels      = [ { 'Volume-based (old version without refinement)' , 'Volume-based (old version with refinement)' 'Surface-based'} , bmethod.labels ];
    bmethod.values      = [ {4, 5}                                                                                                            , bmethod.values ];
  end
  bmethod.val           = {3};
  bmethod.help          = {[ ...
    'Bone processing method using volumes or additional surfaces to extract bone intensities. ' ...
    'Values are normlized for tissue contrast but still depending on image weighting (e.g. T1, T2, PD, EPI) and image protocol parameters (e.g. fat supression) ' ...
    'and harmonization (e.g. via external tools such as COMBAT) is recommendet. '];''};

  % not sure if this is useful
  % SPM/CAT will come with an affine registration that focuses on brain tissues rather than the skull
  % but this should not affect the bone measures alhtough the normalized output is not optimal  
  affreg                = cfg_menu;
  affreg.tag            = 'affreg';
  affreg.name           = 'Apply affine registration';
  affreg.labels         = {'No','Yes'};
  affreg.values         = {0,1};
  affreg.val            = {0};
  affreg.hidden         = expertgui<2; %%%%%%
  affreg.help           = {'Apply bone-focused affine registration.';''};

  rerun                 = cfg_menu;
  rerun.tag            = 'rerun';
  rerun.name           = 'Rerun preprocessing';
  rerun.labels         = {'No','Yes'};
  rerun.values         = {0,1};
  rerun.val            = {0};
  rerun.hidden         = expertgui<2;
  rerun.help           = {'Run processing even if the output already exist.';''};

  reduce                = cfg_menu;
  reduce.tag            = 'reduce';
  reduce.name           = 'Processing resolution';
  reduce.labels         = {'full','half','third','quater'};
  reduce.values         = {1,2,3,4};
  reduce.val            = {2};
  reduce.hidden         = expertgui<1;
  reduce.help           = {
    ['Surface resolution reduction factor as support accurate, robust, and at once fast processing. ' ...
     'Without reduction, surfaces of volumes with about 1 cm resolution typically have 120k vertices in humans. ' ...
     'But for our purpose the reduction to 2 mm with about 30k surfaces is sufficient. ' ...
     'However, reductions for 3 or 4 times result still in quite similare measurements. ']; ''};
 
  mask                  = cfg_files;
  mask.tag              = 'Pmask';
  mask.name             = 'Bone mask';
  mask.help             = {'Bone mask to avoid critical regions not suited for MRI (e.g. facial bones and thin or error prone inferior areas). '};
  mask.filter           = {'image'};
  mask.dir              = fullfile(spm('dir'),'toolbox','boney'); 
  mask.val              = {{fullfile(spm('dir'),'toolbox','boney','boney_KADA_bone-mask.nii')}}; 
  mask.ufilter          = '.*';
  mask.num              = [0 1];
  mask.hidden           = expertgui<1;
  mask.help             = {'Mask problematic regions.';''};
 
  atlas                 = cfg_files;
  atlas.tag             = 'Patlas';
  atlas.name            = 'Bone atlas';
  atlas.help            = {'Bone atlas that describe the major skull bones (frontal, pariatal-left/right, temporal, occipital). '};
  atlas.filter          = {'image'};
  atlas.dir             = fullfile(spm('dir'),'toolbox','boney'); 
  atlas.val             = {{fullfile(spm('dir'),'toolbox','boney','boney_KADA_bone-regions.nii')}}; 
  atlas.ufilter         = '.*';
  atlas.num             = [0 1];
  atlas.hidden          = expertgui<1;

  % main parameter structure passed to cat_plot_boxplot
  opts                  = cfg_exbranch;
  opts.tag              = 'opts';
  opts.name             = 'Options';
  opts.val              = { pmethod , bmethod, affreg, reduce, atlas, mask, nproc, rerun, verb }; 
  opts.help             = {'Specify processing parameters and atlas/mask files.'};



  % == output ==  
  report                = cfg_menu;
  report.tag            = 'report';
  report.name           = 'Write report';
  report.labels         = {'No','Basic','Full'};
  report.values         = {0,2,3};
  report.val            = {3};
  report.hidden         = expertgui<1;
  report.help           = {'Write JPG report file with values only (Basic) or additional volume and available surfaces (Full) into a "report" subdirectory.';''};

  writevol              = cfg_menu;
  writevol.tag          = 'writevol';
  writevol.name         = 'Write volumes';
  writevol.labels       = {'No','Yes'};
  writevol.values       = {0,1};
  writevol.val          = {0};
  writevol.hidden       = expertgui<1;
  writevol.help         = {'Write refined bone segment volumes used for extraction into a "vol" subdirectory';''};
 
  writesurf             = cfg_menu;
  writesurf.tag         = 'writesurf';
  writesurf.name        = 'Write surfaces';
  writesurf.labels      = {'No','Yes'};
  writesurf.values      = {0,1};
  writesurf.val         = {0};
  writesurf.hidden      = expertgui<1;
  writesurf.help        = {'Write (individual) bone surfaces used for extraction (if processed) into a "surf" subdirectory.';
                           'These surfaces cannot be compared directly as they vertices are not aligned.'};

  % output parameter structure
  output                = cfg_exbranch;
  output.tag            = 'output';
  output.name           = 'Output';
  output.val            = { report , writevol , writesurf }; 
  output.help           = {'Specify output parameters. The main results are writen as XML/MAT file into a "report" subdirectory.' ''};



  % == main == 
  segment               = cfg_exbranch;
  segment.tag           = 'segment';
  segment.name          = 'Bone segmentation';
  segment.val           = {files,opts,output}; 
  segment.prog          = @boney_segment;
  segment.vout          = @vout_boney_cfg_segment;
  segment.help          = {
   '' ''};
end
function dep = vout_boney_cfg_segment(job)
%vout_segment. SPM dependency structure for boney_segment. 

  dep            = cfg_dep;
  dep.sname      = 'XML';
  dep.src_output = substruct('.','xml');
  dep.tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});

  if job.output.writevol
    dep(end+1)          = cfg_dep;
    dep(end).sname      = 'Volumes';
    dep(end).src_output = substruct('.','volumes');
    dep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  
  if job.output.writesurf
    dep(end+1)          = cfg_dep;
    dep(end).sname      = 'Surfaces';
    dep(end).src_output = substruct('.','surfaces');
    dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
  end
  
  if job.output.report
    dep(end+1)          = cfg_dep;
    dep(end).sname      = 'Reports';
    dep(end).src_output = substruct('.','reports');
    dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
  end
end




  