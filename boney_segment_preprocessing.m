function P = boney_segment_preprocessing(P,out,ctpm,pmethod,bias,rerun)
%boney_segment_preprocessing. Call SPM12 or CAT12 segmentation.
% 
%  P = boney_segment_preprocessing(P,out,method,bias,rerun)
%
%  P        .. updated structures of processed files 
%  out(:).P .. other input files 
%  ctpm     .. TPM selector (1-default SPM TPM for adults, 2-children)
%  pmethod  .. preprocesing method (1-SPM; 2-CAT; ?-segCT)
%  bias     .. use strong bias correction (i.e., 30 mm)
%  rerun    .. run preprocessing even if all files exist
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________
  
% TODO: 
% * Optimization of the segmentation data storing as one label rather than
%   many segment files. 

  % check for already processed files
  PC = spm_file(P,'ext','.nii'); PC{1} = out(1).P.org; for i=2:numel(PC), PC{i}  = out(i).P.org; end
  if rerun
    PC = P; 
  else
    % extract processed filenames
    Ppc = P; Ppc{1} = out(1).P.cls{1}; for i=2:numel(P), Ppc{i} = out(i).P.cls{1}; end
    Pbc = P; Pbc{1} = out(1).P.bc;     for i=2:numel(P), Pbc{i} = out(i).P.bc{1}; end

    %% have to use CAT developer GUI for rerun function 
    oldexpertgui = cat_get_defaults('extopts.expertgui');
    cat_get_defaults('extopts.expertgui',2);
    rpc = cat_io_rerun(PC,Ppc,0) > 0; 
    rbc = cat_io_rerun(PC,Pbc,0) > 0; % estimate if reprocessing is required
    rsc = rpc | rbc;
    cat_io_cprintf([0 .5 0],'  %d of %d SPM segmentations can be used. \n', numel(rsc) - sum(rsc), numel(rsc));
    if sum(rsc) > 0
        cat_io_cprintf([0.7 .2 0],'  %d of %d files need preprocessing. \n', sum(rsc), numel(rsc));
    end
    
    
    cat_get_defaults('extopts.expertgui',oldexpertgui);

% ################
% * check TPM !!!!
%   we have to load the seg8 of each T1 and check if it was processed with
%   the correct TPM :/
% ################
    
    PC  = PC(rsc); % only keep cases that need preprocessing
  end

  %% processing
  if ~isempty(PC)
    switch ctpm
      case 1, Ptpm = fullfile(spm('dir'),'tpm',sprintf('TPM.nii')); 
      case 2, Ptpm = fullfile(spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym',sprintf('TPM_Age11.5.nii'));
    end
    switch pmethod
      case {1,'spm'}
        matlabbatch = SPM_preprocessing(PC, Ptpm, bias);
      case {2,'cat'} % not working
        expertgui = cat_get_defaults('extopts.expertgui');
        if expertgui < 1, cat12('expert'); end % need expert/developer batch to write class 4 to 6
        matlabbatch = CAT_preprocessing(PC, Ptpm, bias, expertgui); 
      case {3,'CT_seg'}
        if exist(fullfile(spm('dir'),'toolbox','CTseg'),'dir')
          fprintf('CTseg is not working')
          return
        else
          error('CTseg is not available')
        end
      otherwise
        error('unkown preprocessing method')
    end
  
    % data optimization 
% #####################
% The idea was to combine all 6 classes in one label map to save space but 
% also time in loading these images from disk (factor 6!). 
% However, this would need further additional checks and case handling.
% #####################

    % run SPM batch
    warning off; 
    spm_jobman('run',matlabbatch);
    warning on; 

  end


end
% subfunction with the main CAT matlabbatch
function matlabbatch = SPM_preprocessing(P,Ptmp,bias)
%SPM_preprocessing. Create SPM12 segmentation matlabbatch.
  matlabbatch{1}.spm.spatial.preproc.channel.vols        = P;
  matlabbatch{1}.spm.spatial.preproc.channel.biasreg     = 0.001;
  matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm    = 60 - 30*bias;
  matlabbatch{1}.spm.spatial.preproc.channel.write       = [0 1]; % write bias corrected image
  % use default classes (a skull-stripped case is irrelevant here ;-)
  ngaus = [1 1 2 3 4 2]; 
  for ti = 1:6
    matlabbatch{1}.spm.spatial.preproc.tissue(ti).tpm    = {sprintf('%s,%d',Ptmp,ti)};
    matlabbatch{1}.spm.spatial.preproc.tissue(ti).ngaus  = ngaus(ti);
    matlabbatch{1}.spm.spatial.preproc.tissue(ti).native = [ti<6 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(ti).warped = [0 0];
  end
  matlabbatch{1}.spm.spatial.preproc.warp.mrf            = 1;
  matlabbatch{1}.spm.spatial.preproc.warp.cleanup        = 1;
  matlabbatch{1}.spm.spatial.preproc.warp.reg            = [0 0.001 0.5 0.05 0.2];
  % subj was more robust compair to the default (=mri) 
  matlabbatch{1}.spm.spatial.preproc.warp.affreg         = 'subj';
  matlabbatch{1}.spm.spatial.preproc.warp.fwhm           = 0;
  % surprisingly samp=5 was more robust than 3 which failed in good bone segmenation in many UKB cases
  matlabbatch{1}.spm.spatial.preproc.warp.samp           = 5; % default is 3 
  matlabbatch{1}.spm.spatial.preproc.warp.write          = [0 0];
  matlabbatch{1}.spm.spatial.preproc.warp.vox            = NaN;
  matlabbatch{1}.spm.spatial.preproc.warp.bb             = [NaN NaN NaN; NaN NaN NaN];
end
% subfunction with the main CAT matlabbatch
function matlabbatch = CAT_preprocessing(P, Ptmp, bias, expertgui)
%CAT_preprocessing. Create CAT12 segmentation matlabbatch.
  matlabbatch{1}.spm.tools.cat.estwrite.data                              = P;
  matlabbatch{1}.spm.tools.cat.estwrite.data_wmh                          = {''};
  matlabbatch{1}.spm.tools.cat.estwrite.nproc                             = 0; % if parallel, then maybe the whole processing
  matlabbatch{1}.spm.tools.cat.estwrite.useprior                          = '';
  matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm                          = {Ptmp};
  matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg                       = 'subj';
  if expertgui == 2
    matlabbatch{1}.spm.tools.cat.estwrite.opts.ngaus                      = [1 1 2 3 4 2];
    matlabbatch{1}.spm.tools.cat.estwrite.opts.warpreg                    = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.tools.cat.estwrite.opts.redspmres                  = 0;
  end
  matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr                      = 0.5 + 0.5*bias;
  matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr                       = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.restypes.optimal = [1 0.3];
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.setCOM       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP          = 1070;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.affmod       = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr        = .5; % faster? - no, better robust 
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.spm_kamap    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASstr       = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASmyostr    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr      = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr   = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.BVCstr       = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHC         = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.SLC          = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.mrf          = 1;
  if expertgui == 2
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.T1           = ...
      {fullfile(spm('dir'), 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym/T1.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.brainmask    = ...
      {fullfile(spm('dir'), 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym/brainmask.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.cat12atlas   = ...
      {fullfile(spm('dir'), 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym/cat.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.darteltpm    = ...
      {fullfile(spm('dir'), 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym/Template_1_Dartel.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shootingtpm  = ...
      {fullfile(spm('dir'), 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym/Template_0_GS.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regstr       = 0.5;
  else
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regmethod.shooting.shootingtpm = ...
      {fullfile(spm('dir'), 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym/Template_0_GS.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regmethod.shooting.regstr = 0.5;
  end
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.bb           = 12;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.vox          = 1.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtres            = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtmethod         = 'pbt2x';
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.SRP               = 30;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.reduce_mesh       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.vdist             = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.scale_cortex      = 0.7;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.add_parahipp      = 0.1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.close_parahipp    = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.experimental        = 0; 
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.new_release         = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy                = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors        = 1; %
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb                = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print               = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno                = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.surface                    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures              = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics                       = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40                                   = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra                                    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers                                  = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus                                 = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamic_nuclei                          = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.suit                                     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr                                     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal3                                     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori                                     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy3                                 = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.julichbrain                              = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Tian_Subcortex_S4_7T                     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_100Parcels_17Networks_order = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_200Parcels_17Networks_order = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_400Parcels_17Networks_order = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_600Parcels_17Networks_order = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas                                 = {''};
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native      = 1; % GM
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod         = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native      = 1; % WM
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod         = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native     = 1; % CSF
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod        = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.pp.native      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.pp.warped      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.pp.dartel      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod        = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod         = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native    = 1; % head
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod       = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native   = 0;
  if expertgui == 2
    matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.dartel = 0;
  end
  matlabbatch{1}.spm.tools.cat.estwrite.output.label.native   = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped   = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel   = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native    = 1; % image (m*)
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped    = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel    = 0; 
  matlabbatch{1}.spm.tools.cat.estwrite.output.las.native     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel     = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.warps          = [1 0]; % normalization ???
  matlabbatch{1}.spm.tools.cat.estwrite.output.rmat           = 0;
end
