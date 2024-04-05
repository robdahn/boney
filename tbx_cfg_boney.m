function boney = tbx_cfg_boney
% Configuration file for boney toolbox batches.
% _________________________________________________________________________
% ###### Please write more about this ...
% _________________________________________________________________________
% Robert Dahnke 2023
%
%#ok<*AGROW,*INUSD>

  global boned %#ok<GVMIS> 

  if isfield(boned,'expertgui')
    expertgui = boned.expertgui; 
  else
    expertgui = 2; 
  end
  
  addpath(fileparts(which(mfilename)));
  
  
  % try to estimate the number of processor's cores
  try
    numcores = cat_get_defaults('extopts.nproc');
    % use only half of the cores for Windows due to inadequate memory management 
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
      'In order to use multi-threading jobs with multiple subjects processing can be split into separate processes that run in the background. ' ...
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
  files.help          = {'Select images that should be processed or select the bias-corrected processed images "m*.nii".'};
  files.filter        = 'image';
  files.ufilter       = '.*';
  files.num           = [1 Inf];
   
  verb                = cfg_menu;
  verb.tag            = 'verb';
  verb.name           = 'Verbose processing (expert)';
  verb.labels         = {'No','Yes'};
  verb.values         = {0,1};
  if expertgui 
    verb.labels       = [ verb.labels, {'Yes - Details'} ];
    verb.values       = [ verb.values, {2} ];
    verb.help         = {'Verbose processing. Currently with details to support more information in testing. '};
  end
  verb.val            = {1}; % RD202403: avoid details {max(1,expertgui)}; 
  verb.hidden         = expertgui<1; 
  
  
  %  batches
  %  ------------------------------------------------------------------------
  segment = boney_cfg_segment(files,nproc,expertgui,verb);
  xml2csv = conf_io_xml2csv(expertgui); 
  
  %  main
  %  ------------------------------------------------------------------------
  boney               = cfg_choice;
  boney.name          = 'Boney';
  boney.tag           = 'boney';
  boney.values        = {segment;xml2csv};
return


%  subfunctions
%  ------------------------------------------------------------------------
function segment = boney_cfg_segment(files,nproc,expertgui,verb)
%boney_cfg_segment. Bone segmentation and value extraction.

  % == opts == 
  % development option - start with SPM, will add CAT which is maybe more
  % robust in some cases, maybe use different presettings, e.g. TPMs for
  % children later 
  pmethod               = cfg_menu;
  pmethod.tag           = 'pmethod';
  pmethod.name          = 'Preprocessing method';
  pmethod.labels        = {'SPM','CAT'}; % ###### add CTseg later in case the toolbox is installed.  
  pmethod.values        = {1,2};
  pmethod.val           = {1};
  pmethod.help          = {['Preprocessing method to segment different tissue classes in the given image. ' ... 
    'For SPM, the bias-corrected input "m*.nii" is used (intensity normalized for BG-WM), '...
    'whereas for CAT12 the "m*.nii" is a denoised, bias-corrected, and step-wise intensity normalized image (with BG=0,CSF=1/3,GM=2/3,WM=1).' ...
    ' Therefore, the intensity-based values are not comparable, whereas the thickness-based values should be similar.'];''};

  % main method of this toolbox
  bmethod               = cfg_menu;
  bmethod.tag           = 'bmethod';
  bmethod.name          = 'Bone processing method';
  bmethod.labels        = {'SPM mat-file', 'Volume-based', 'Surface-based'};
  bmethod.values        = {0, 1, 2};
  %########################## not prepared yet
  %{
  if expertgui 
    bmethod.labels      = [ bmethod.labels , { 'Volume-based (old version)' } ];
    bmethod.values      = [ bmethod.values , { 3 } ];
  end
  %}
  bmethod.val           = {2};
  bmethod.help          = {[ ...
    'Bone processing method using volumes or additional surfaces to extract bone intensities. ' ...
    'Values are normalized for tissue contrast but are still depending on image weighting (e.g., T1, T2, PD, EPI) and image protocol parameters (e.g., fat supression). ' ...
    'Harmonization (e.g., via external tools such as COMBAT) is recommended. ']; ''};

  classic               = cfg_menu;
  classic.tag           = 'classic';
  classic.name          = 'Estimate also fast classic/prototype measures (developer)';
  classic.labels        = {'No','Yes'};
  classic.values        = {0,1};
  classic.val           = {expertgui==2};
  classic.hidden        = expertgui<2;
  classic.help          = {'Estimate also fast classic measures (prototype).';''};

  % RD202309: not sure if this is useful
  % * SPM/CAT will come with an affine registration that focuses on brain 
  %   tissues rather than the skull but this should not affect the bone 
  %   measures as much, although the normalized output is not optimal  
  % * registration could be interesting with an additional non-linear Shooting
  %   that would need extra optimized feature classes, e.g. separate classes 
  %   for lobes or other head features
  % * maybe a surface-based registration would work better here
  %   >> would need manual segmentation to test this
  % * However, it is still unclear if this would support some better output
  %   Maybe for machine or deep learning. 
  affreg                = cfg_menu;
  affreg.tag            = 'affreg';
  affreg.name           = 'Apply affine registration (developer)';
  affreg.labels         = {'No','Yes'};
  affreg.values         = {0,1};
  affreg.val            = {0};
  affreg.hidden         = expertgui<3; % ########### not implemented yet
  affreg.help           = {'Apply bone-focused affine registration.';''};

  refine                = cfg_menu;
  refine.tag            = 'refine';
  refine.name           = 'Refine preprocessing (expert)';
  refine.labels         = {'No','Yes','Yes - enhanced'};
  refine.values         = {0,1,2};
  refine.val            = {2};
  refine.hidden         = expertgui<1;
  refine.help           = {[ ...
    'Without fat supression, the bone marrow can have high intensities that can be mislabed as head. ' ...
    'This can be seen as large local underestimations of bone thickness and bone intensity - ' ...
    'typically well visible on the skull surface. ' ...
    'Morphological operations were used to close such holes and obtain a more complete skull segment. '];'';
   ['The enhanced version made a precorrection of the background and brainmask to estimate a distance-based refinement based on atlas regions. ' ...
    'The assumption is that the head and skull have a 50:50 relationship in thinner regions of the skull, which would be covered by a hat, to correct severe outliers. '];'';
   };
  
  % As this is seen to be stable - reprocessing is generally not required.
  % To support faster reprocessing of the bone measures, it would be nice
  % to store the segmentation files - maybe compressed? 
  % Anyhow, this is an expert feature. 
  prerun                = cfg_menu;
  prerun.tag            = 'prerun';
  prerun.name           = 'Rerun SPM/CAT preprocessing (expert)';
  prerun.labels         = {'No','Yes'};
  prerun.values         = {0,1};
  prerun.val            = {0};
  prerun.hidden         = expertgui<1;
  prerun.help           = { ...
    ['Run SPM/CAT processing even if the output already exists. ' ...
     'As preprocessing is considered to be stable, reprocessing is generally not required. '];''};

  rerun                 = cfg_menu;
  rerun.tag             = 'rerun';
  rerun.name            = 'Rerun bone processing (expert)';
  rerun.labels          = {'No','Yes'};
  rerun.values          = {0,1};
  rerun.val             = {1};
  rerun.hidden          = expertgui<1;
  rerun.help            = {'Run bone processing even if the output already exists.';''};

  % General resolution limit:   
  % Would be good to use this even before SPM/CAT to stabilize and speed up  
  % the whole preprocessing. However, then the output images would have 
  % also a different size. 
  reslim                = cfg_menu;
  reslim.tag            = 'reslim';
  reslim.name           = 'Resolution limit (expert)';
  reslim.labels         = {'0.5 mm','1.0 mm', '1.5 mm'};
  reslim.values         = {0.5,1,1.5};
  reslim.val            = {1.5};
  reslim.hidden         = expertgui<1;
  reslim.help           = ... 
   {['Limit processing resolution for the whole pipeline to stabilize and speed-up the processing. ' ...
     'The limit gives the required resolution after reduction, i.e., an image of .75 would be reduced ' ...
     'to 1.5 mm, whereas 0.8 would not be reduced. '];''};
  
  reduce                = cfg_menu;
  reduce.tag            = 'reduce';
  reduce.name           = 'Surface processing resolution (expert)';
  reduce.labels         = {'Full','Half','Third','Quarter'};
  reduce.values         = {1,2,3,4};
  reduce.val            = {2};
  reduce.hidden         = expertgui<1;
  reduce.help           = {
    ['The reduction factor of surface resolution supports accurate, robust, and fast processing. ' ...
     'Without reduction, surfaces of volumes with about 1 cm resolution typically have 120k vertices in humans. ' ...
     'For our purpose, the reduction to 2 mm with about 30k surfaces is sufficient. ' ...
     'However, reductions for 3 or 4 times still result in quite similar measurements. ']; ''};
 
  % The idea was to avoid problematic regions in the global case. 
  % However, this cannot be applied to the atlas regions in general.
  mask                  = cfg_files;
  mask.tag              = 'Pmask';
  mask.name             = 'Bone mask (expert)';
  mask.help             = {'Bone mask to avoid critical regions not suited for MRI analysis (e.g., facial bones and thin or error prone areas). '};
  mask.filter           = {'image'};
  mask.dir              = fullfile(spm('dir'),'toolbox','boney'); 
  mask.val              = {{fullfile(spm('dir'),'toolbox','boney','boney_KADA_bone-mask.nii')}}; 
  mask.ufilter          = '.*';
  mask.num              = [0 1];
  mask.hidden           = expertgui<1; 
  mask.help             = {'Mask problematic regions for global but not regional measurements. Deselect a file to avoid masking. ';''};
 
  % Would it be useful to apply cortical atlases here? 
  % - Difficult, as only regions to the skull would be mapped and further
  %   differences for SPM/CAT would come by the registration.
  % Would a finer regional resolution be useful?
  % - Maybe, if the regions have some useful definition.
  % - Could be also some regularly defined subregions by some kind of parcelation. 
  % - This (especially finer regions) would perhaps require some non-linear registration. 
  % - However, here we would come to VBM, ML or DL, which are also
  %   connected to the general registration question.
  % Hence, the altas is an expert feature
  atlas                 = cfg_files;
  atlas.tag             = 'Patlas';
  atlas.name            = 'Bone atlas (expert)';
  atlas.help            = { ...
    ['Bone atlas that describes the major skull bones (e.g., frontal, occipital, parietal-left/right, temporal-left/right). ' ...
     'The atlas is mapped by the affine registration from the SPM/CAT preprocessing and linearly extended. ' ...
     'The overlap between sutures and atlas boundaries is adequate, but could be improved in future. ' ...
     'As we only have 6 large major regions, regional extraction seems to be sufficient so far. '];''};
  atlas.filter          = {'image'};
  atlas.dir             = fullfile(spm('dir'),'toolbox','boney'); 
  atlas.val             = {{fullfile(spm('dir'),'toolbox','boney','boney_KADA_bone-regions.nii')}}; 
  atlas.ufilter         = '.*';
  atlas.num             = [0 1];
  atlas.hidden          = expertgui<1;

  % I want to keep the menu simple and therefore avoid a more general file selector
  ctpm                  = cfg_menu;
  ctpm.tag              = 'ctpm';
  ctpm.name             = 'TPM';
  ctpm.labels           = {'Adults','Childen'}; % some additional small-child template would be nice 
  ctpm.values           = {1,2};
  ctpm.val              = {1};
  ctpm.help             = {'Select TPM (Tissue Probility Map) for preprocessing. ';''};

  
  % main parameter structure passed to cat_plot_boxplot
  opts                  = cfg_exbranch;
  opts.tag              = 'opts';
  opts.name             = 'Options';
  opts.val              = { ctpm, pmethod, bmethod, classic, reslim, affreg, ...
    refine, reduce, atlas, mask, nproc, prerun, rerun, verb }; 
  opts.help             = {'Specify processing parameters and atlas/mask files. '};



  % == output ==  
  % this is probably the only flag that a user can set without risk that
  % could result in saving processing speed 
  report                = cfg_menu;
  report.tag            = 'report';
  report.name           = 'Write report';
  report.labels         = {'No','Basic','Full'};
  report.values         = {0,2,3};
  report.val            = {3};
  %report.hidden         = expertgui<1;
  report.help           = {['Write JPG report file with values only (Basic) or ' ...
    'additional volume and available surfaces (Full) into a "report" subdirectory.'];''};

  resdir                = cfg_entry;
  resdir.tag            = 'resdir';
  resdir.name           = '(Relative) folder';
  resdir.strtype        = 's';
  resdir.num            = [0 Inf];
  resdir.val            = {'../derivatives/boney'};
  resdir.hidden         = expertgui<1;
  resdir.help           = {
   ['Use relative directory structure to store data, even if it does not conform to BIDS.  ' ...
    'This alternative definion, based on the depth of the file and controlled here by the repetition of "../", ' ...
    'is keeping the structure of subdirectories more robust in case of a regular, but non-BIDS structure without default directory, ' ...
    'naming "sub-##/ses-##/anat" and similar filenames, e.g. for "../../derivatives/CAT##.#_#" and the following files:'];
    '   ../group-01/sub-01/t1w.nii';
    '   ../group-01/sub-02/t1w.nii';
    'it results in:'; 
    '   ../derivatives/CAT##.#_#/group-01/sub-01/t1w.nii';
    '   ../derivatives/CAT##.#_#/group-01/sub-02/t1w.nii';
    'rather than:';
    '   ../derivatives/CAT##.#_#/t1w.nii';
    '   ../derivatives/CAT##.#_#/t1w.nii';
    'where the relative BIDS folder would also cause conflicts by overwriting results.';
    '';
    };

  % alternative naming: delete temporary segmetation files
  % - not sure if data can be compressed in one image - maybe as shell
  %   model with PVE to have just one volume (should be fine enough)
  % - but would also work mat file (that should be compressed)
  writeseg              = cfg_menu;
  writeseg.tag          = 'writeseg';
  writeseg.name         = 'Keep SPM/CAT segmentation (expert)';
  writeseg.labels       = {'No','Yes (original NIFTIs)','Yes (compressed as .mat)'};
  writeseg.values       = {0,1,2};
  writeseg.val          = {expertgui};  
  writeseg.hidden       = expertgui<1;
  writeseg.help         = { ...
   ['Do not delete the native SPM/CAT segmentation files c# mri/p# to support faster processing. ' ....
    'The compressed option saved the optimized volumes in a .mat file to save disk space and loading speed. '];''};

  % - bone & skull ?
  % - settings for native/affine - affine could be used for ML/DL 
  writevol              = cfg_menu;
  writevol.tag          = 'writevol';
  writevol.name         = 'Write bone (expert)';
  writevol.labels       = {'No','Native','Warped','Affine','Native + Affine + Warped'};
  writevol.values       = {0,1,2,3,4};
  writevol.val          = {0};
  writevol.hidden       = expertgui<1;
  writevol.help         = {'Write refined bone segment volumes used for extraction into a "vol" subdirectory.';''};
 
  % - here we have the general registration issue 
  writesurf             = cfg_menu;
  writesurf.tag         = 'writesurf';
  writesurf.name        = 'Write surfaces (expert)';
  writesurf.labels      = {'No','Yes'};
  writesurf.values      = {0,1};
  writesurf.val         = {0};
  writesurf.hidden      = expertgui<1;
  writesurf.help        = {'Write (individual) bone surfaces used for extraction (if processed) into a "surf" subdirectory.';
                           'These surfaces cannot be compared directly as the vertices are not aligned.'};

  % output parameter structure
  output                = cfg_exbranch;
  output.tag            = 'output';
  output.name           = 'Output';
  output.val            = { resdir , report , writeseg, writevol , writesurf }; 
  output.help           = {'Specify output parameters. The main results are written as XML/MAT file into a "report" subdirectory.' ''};



  % == main == 
  segment               = cfg_exbranch;
  segment.tag           = 'segment';
  segment.name          = 'Bone segmentation';
  segment.val           = {files,opts,output}; 
  segment.prog          = @boney_segment;
  segment.vout          = @vout_boney_cfg_segment;
  segment.help          = {
   '' ''};
return
function xml2csv = conf_io_xml2csv(expertgui)
% -------------------------------------------------------------------------
% Read structures/XML-files and export and transform it to a table/CSV-file
% 
% RD202304
% -------------------------------------------------------------------------
  
  % n-files, e.g. XML for direct extraction or nii/gii as selector
  files               = cfg_files;
  files.num           = [1 Inf];
  files.tag           = 'files';
  files.name          = 'XML files';
  files.filter        = 'any';
  files.ufilter       = '^boney_.*\.xml$';
  files.val           = {{''}};
  files.help          = {'Select XML files of one type (e.g., "boney_", "cat_", "catROI_" or "catROIs"). '};

  outdir            = cfg_files;
  outdir.tag        = 'outdir';
  outdir.name       = 'Output directory';
  outdir.filter     = 'dir';
  outdir.ufilter    = '.*';
  outdir.num        = [0 1];
  outdir.help       = {'Select a directory where files are written. Use current directory if empty.'};
  outdir.val{1}     = {''};


  % filename
  fname             = cfg_entry;
  fname.tag         = 'fname';
  fname.name        = 'Filename';
  fname.strtype     = 's';
  fname.num         = [1 inf];
  if expertgui
    fname.val       = {'Boney_xml_REPORT_DATE.csv'}; 
  else
    fname.val       = {'Boney_xml_DATE.csv'}; 
  end
  fname.help        = {'CSV filename.' };
  
  % expert output options
  % fieldnames
  fieldnames            = cfg_entry;
  fieldnames.tag        = 'fieldnames';
  fieldnames.name       = 'Included fieldnames';
  fieldnames.strtype    = 's+';
  fieldnames.val        = {{' '}};
  fieldnames.num        = [0 inf];
  fieldnames.hidden     = expertgui<1; 
  fieldnames.help       = {
    ['Define keywords or complete fields to limit the extraction of the fields (empty = include all). ' ...
     'Only fields that include these strings will be used. ' ...
     'In case of catROI-files you can limit the extraction to specific atlases, regions, or tissues. ' ...
     'In case of cat-files you can limit the extraction to specific parameters ("opts" or "extopts") ' ...
     'or QC ratings ("qualityratings"). ' ...
     'The keyword DATE is replace by the data string "YYYYmmDD-HHMMSS".']
     ''
    };
  if expertgui
    fieldnames.help = [fieldnames.help(1:end-1);{'The keyword REPORT is replace by the selected export level.';''}]; 
  end

  % avoidfields
  avoidfields            = cfg_entry;
  avoidfields.tag        = 'avoidfields';
  avoidfields.name       = 'Excluded fieldnames';
  avoidfields.strtype    = 's+';
  avoidfields.val        = {{''}};
  avoidfields.num        = [0 inf];
  avoidfields.hidden     = expertgui<1; 
  avoidfields.help       = {
     'Define keywords or complete fields that should be avoided/excluded (empty = exclude none). Fields that include these strings will not be exported, even if they were included before. '
     'In case of catROI-files you can limit the extraction to specific atlas, regions, or tissues. '
     'In case of cat-files you can limit the extraction to specific parameters or measures. '
     ''
    };


  report           = cfg_menu;
  report.tag       = 'report';
  report.name      = 'Boney XML export field sets';
  report.labels    = {'Boney default', 'Boney details', 'Boney expert', ...
                      ...'Only processing parameters','No processing parameters', ...
                      'All fields'};
  report.values    = {'boney_default', 'boney_details', 'boney_expert', ...
                      ...'paraonly', 'nopara', ...
                      'all'
                      };
  report.val       = {'boney_default'}; 
  report.help      = {'Predefined sets of boney XML values in case of "boney_" processing XML files (no effect in other XMLs). '};

  
  xml2csv           = cfg_exbranch;
  xml2csv.tag       = 'xml2csv';
  xml2csv.name      = 'XML2CSV';
  xml2csv.val       = {files outdir fname fieldnames avoidfields report};
  xml2csv.prog      = @boney_xml2csv;
  xml2csv.vout      = @vout_boney_cfg_xml2csv; 
  xml2csv.help      = {
    'Export XML files (e.g., the boney preprocessing results) as a CSV table. '
    };

return
function dep = vout_boney_cfg_segment(job)
%vout_boney_cfg_segment. SPM dependency structure for boney_segment. 

  % the XML is the most relevant one
  dep            = cfg_dep;
  dep.sname      = 'XML';
  dep.src_output = substruct('.','xml','()',{':'});
  dep.tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}}); % any cfg_files


  % also the report is generally available 
  if job.output.report
    dep(end+1)          = cfg_dep;
    dep(end).sname      = 'Reports';
    dep(end).src_output = substruct('.','reports');
    dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
  end

  if job.output.writeseg % only native 
    for ci = 1:5
      dep(end+1)          = cfg_dep;
      dep(end).sname      = sprintf('cls %d',ci);
      dep(end).src_output = substruct('.','cls','{}',{'1'},'{}',{':'});
      dep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
  end
  
  % bone volume that includes both cortex and marrow
  if job.output.writevol == 1 || job.output.writevol == 3
    dep(end+1)          = cfg_dep;
    dep(end).sname      = 'Bone(native)';
    dep(end).src_output = substruct('.','bone_native','{}',{':'});
    dep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if job.output.writevol == 2 || job.output.writevol == 3
    dep(end+1)          = cfg_dep;
    dep(end).sname      = 'Bone(affine)';
    dep(end).src_output = substruct('.','rbone_affine','{}',{':'});
    dep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  

  % there are multiple surfaces, but as no registration was done they cannot
  % be used so far (e.g., to create an average etc.)
  if job.output.writesurf
    dep(end+1)          = cfg_dep;
    dep(end).sname      = 'Surfaces';
    dep(end).src_output = substruct('.','bcentral');
    dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});

    dep(end+1)          = cfg_dep;
    dep(end).sname      = 'Surfaces';
    dep(end).src_output = substruct('.','bcentral');
    dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});

    dep(end+1)          = cfg_dep;
    dep(end).sname      = 'Surfaces';
    dep(end).src_output = substruct('.','bmarrow');
    dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});

    dep(end+1)          = cfg_dep;
    dep(end).sname      = 'Surfaces';
    dep(end).src_output = substruct('.','bthickness');
    dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});

    dep(end+1)          = cfg_dep;
    dep(end).sname      = 'Surfaces';
    dep(end).src_output = substruct('.','hthickness');
    dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
  end
  
return
function dep = vout_boney_cfg_xml2csv(job)
%vout_boney_cfg_xml2csv. SPM dependency structure for boney_xml2csv. 

  dep            = cfg_dep;
  dep.sname      = 'CSV';
  dep.src_output = substruct('.','csv','()',{':'});
  dep.tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}}); % any cfg_files
return



  
