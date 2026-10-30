function varargout = boney_segment_preprocessing(P,out,job,outbc,Ptpmfb)
%boney_segment_preprocessing. Call SPM12 or CAT12 segmentation.
%
%  boney_segment_preprocessing(P,out,job[,outbc,Ptpmfb])
%
%  P        .. images that have to be processed
%  out(:).P .. other input files
%  job      .. main job structure with
%   .opts.ctpm    .. TPM selector (0-automatic, 1-default SPM TPM for adults,
%                    2-children, 3-BlaiottaTPM, 4-optimized TPM7blr)
%   .opts.pmethod .. preprocesing method (1-SPM; 2-CAT; ?-segCT)
%   .opts.bias    .. use strong bias correction (i.e., 40 mm)
%   .opts.spm     .. SPM parameters (see boney_segment)
%  outbc    .. subjects with SPM segmentation but without bias corrected
%               image that is written by the existing SPM parameters
%  Ptpmfb   .. TPM of the automatic fallback (see boney_segment) that is
%               used instead of the automatic TPM selection
%
% SPM is called with one batch for all subjects of a TPM. If this batch 
% fails (e.g. due to a corrupt image), the subjects are processed 
% separately so that only the problematic subject has no segmentation 
% (and gets an error in the main processing).
%
% Further calls to get default values:
%  Ptpms = boney_segment_preprocessing('tpms');        % TPMs of ctpm 1-4
%  ngaus = boney_segment_preprocessing('ngaus',Ptpm);  % default Gaussians
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

  % default values
  if ischar(P)
    switch P
      case 'tpms'
        % the 4th TPM is optional (only some SPM installations) and not
        % used by the automatic selection
        varargout{1} = {
          fullfile(spm('dir'),'tpm','TPM.nii')
          fullfile(spm_file(which('cat12'),'fpath'),'templates_MNI152NLin2009cAsym','TPM_Age11.5.nii')
          fullfile(spm('dir'),'tpm','BlaiottaTPM.nii')
          fullfile(spm('dir'),'tpm','mni0R1p5_TPM7blr.nii')};
      case 'ngaus'
        % 7-class TPMs with the BlaiottaTPM class order (GM, WM, CSF,
        % soft tissue, fat, skull, background) and 6-class SPM TPMs
        if numel(spm_vol(spm_file(out,'number',''))) == 7
          varargout{1} = [1 1 1 3 1 2 3];
        else
          varargout{1} = [1 1 2 3 4 2];
        end
    end
    return
  end

  if ~exist('outbc','var'), outbc = []; end

  % write the bias corrected image if only this one is missing
  for si = 1:numel(outbc)
    writeBiasCorrected(outbc(si));
  end

  if isempty(P), return; end

  %% processing
  Ptpms = boney_segment_preprocessing('tpms');
  switch job.opts.pmethod
    case {1,'spm'}
      % unzip gz files temporary (the original file is kept)
      PCS = P; delunzipped = false(size(P));
      for pi = 1:numel(P)
        [pp,ff,ee] = spm_fileparts(P{pi});
        if strcmp(ee,'.gz')
          PCS{pi} = fullfile(pp,ff);
          if ~exist(PCS{pi},'file')
            delunzipped(pi) = true;
            gunzip(P{pi});
          end
        end
      end

      % TPM for each subject
      Ptpm = cell(size(P));
      for pi = 1:numel(P)
        if exist('Ptpmfb','var') && ~isempty(Ptpmfb)
          Ptpm{pi} = Ptpmfb; % the selection file is updated by boney_segment
        elseif ~isempty(job.opts.spm.tpm) || job.opts.ctpm > 0
          if ~isempty(job.opts.spm.tpm)
            Ptpm{pi} = job.opts.spm.tpm;
          else
            Ptpm{pi} = Ptpms{job.opts.ctpm};
          end
          if ~exist(Ptpm{pi},'file')
            error('boney:preprocessing:missingTPM','TPM "%s" is not available.',Ptpm{pi});
          end
          % remove the selection of a previous automatic run
          if exist(out(pi).P.tpmsel,'file'), delete(out(pi).P.tpmsel); end
        else
          % automatic TPM selection (a corrupt image fails already here)
          try
            [ctpm, tpmsel] = boney_segment_selectTPM(PCS{pi}, job.opts.verb);
          catch e
            reportFailed(PCS(pi),e);
            continue
          end
          tpmsel.fallback = 0;
          tpmsel.tried    = ctpm;
          Ptpm{pi} = Ptpms{ctpm};
          save(out(pi).P.tpmsel,'-struct','tpmsel');
        end
      end

      % run SPM batch (one batch for each TPM)
      valid = ~cellfun('isempty',Ptpm);
      [uTPM,~,tpmi] = unique(Ptpm(valid)); PCSv = PCS(valid);
      for ti = 1:numel(uTPM)
        runSPMbatch(PCSv(tpmi==ti), uTPM{ti}, job.opts);
      end

      for pi = 1:numel(P)
        if delunzipped(pi), delete( PCS{pi} ); end
      end

    case {2,'cat'}
      if job.opts.ctpm == 0 || job.opts.ctpm >= 3
        Ptpm = Ptpms{1}; % CAT requires a 6 class TPM
      else
        Ptpm = Ptpms{job.opts.ctpm};
      end
      expertgui = cat_get_defaults('extopts.expertgui');
      guinames  = {'default','expert','developer'};
      if expertgui < 1, oldexpertgui = expertgui; cat12('expert'); end % need expert/developer batch to write class 4 to 6
      matlabbatch = CAT_preprocessing(P, Ptpm, job.opts.bias, expertgui);
      if exist('oldexpertgui','var'), cat12(guinames{oldexpertgui+1}); end
      warning off;
      spm_jobman('run',matlabbatch);
      warning on;

    case {3,'CT_seg'}
      if exist(fullfile(spm('dir'),'toolbox','CTseg'),'dir')
        error('boney:preprocessing:CTseg','CTseg preprocessing is not supported yet.')
      else
        error('boney:preprocessing:CTseg','CTseg is not available.')
      end
    otherwise
      error('boney:preprocessing:method','Unknown preprocessing method.')
  end

end
%=======================================================================
function runSPMbatch(P,Ptpm,opts)
%runSPMbatch. Run one SPM batch for all images. If it fails (e.g. due to 
% a corrupt image), the images are processed separately so that only the
% problematic images remain without segmentation. Images that were 
% already segmented by the failed batch are not processed again.
  matlabbatch = SPM_preprocessing(P, Ptpm, opts);
  t0 = now; %#ok<TNOW1>
  try
    warning off;
    spm_jobman('run',matlabbatch);
    warning on;
  catch e
    warning on;
    if numel(P) == 1, reportFailed(P,e); return; end
    cat_io_cprintf('warn',sprintf(['  SPM segmentation of %d images failed (%s).\n' ...
      '  Process the images separately.\n'], numel(P), e.message));
    for pi = 1:numel(P)
      if isNewFile(spm_file(P{pi},'suffix','_seg8','ext','.mat'),t0) && ...
         isNewFile(spm_file(P{pi},'prefix','c1'),t0)
        continue
      end
      try
        matlabbatch = SPM_preprocessing(P(pi), Ptpm, opts);
        warning off;
        spm_jobman('run',matlabbatch);
        warning on;
      catch e
        warning on;
        reportFailed(P(pi),e)
      end
    end
  end
end
%=======================================================================
function reportFailed(P,e)
%reportFailed. Report a failed SPM segmentation (or TPM selection) but 
% proceed with the other images (the main processing reports an error for
% this image as its segmentation is missing).
  cat_io_cprintf('err',sprintf('  SPM segmentation failed for "%s": %s\n', P{1}, e.message));
end
%=======================================================================
function isnew = isNewFile(P,t0)
%isNewFile. True if the file exists and was written after t0.
  d     = dir(P);
  isnew = ~isempty(d) && d(1).datenum >= t0;
end
%=======================================================================
function writeBiasCorrected(out)
%writeBiasCorrected. Write the bias corrected image m* based on the bias
% field of an existing SPM segmentation (e.g. of a default SPM run that
% does not write the m* image) without changing the segmentation.
% The bias field is reconstructed as in spm_preproc_write8.

  s8 = load(out.P.seg8,'Tbias','image');
  [pp,ff,ee] = spm_fileparts(out.P.org);
  if strcmp(ee,'.gz')
    Porg = fullfile(pp,ff);
    delunzipped = ~exist(Porg,'file');
    if delunzipped, gunzip(out.P.org); end
  else
    Porg = out.P.org; delunzipped = false;
  end
  V    = spm_vol(Porg);
  d    = V.dim(1:3);
  T    = s8.Tbias{1};
  d3   = size(T); d3(end+1:3) = 1;
  B1   = spm_dctmtx(d(1),d3(1));
  B2   = spm_dctmtx(d(2),d3(2));
  B3   = spm_dctmtx(d(3),d3(3));
  Ym   = zeros(d,'single');
  for z = 1:d(3)
    t1  = reshape(reshape(T, d3(1)*d3(2), d3(3)) * B3(z,:)', d3(1), d3(2));
    bf  = exp(B1 * t1 * B2');
    Ym(:,:,z) = single(spm_slice_vol(V, spm_matrix([0 0 z]), d(1:2), 0)) .* bf;
  end
  Vm         = rmfield(V,{'private'});
  Vm.fname   = out.P.bc;
  Vm.dt      = [spm_type('float32') spm_platform('bigend')];
  Vm.pinfo   = [1;0;0];
  Vm.descrip = 'Bias corrected (reconstructed by boney)';
  spm_write_vol(Vm,Ym);
  if delunzipped, delete(Porg); end
end
%=======================================================================
function matlabbatch = SPM_preprocessing(P,Ptpm,opts)
%SPM_preprocessing. Create SPM12 segmentation matlabbatch.
  if isempty(opts.spm.biasfwhm)
    biasfwhm = 60 - 20*opts.bias; % strong bias correction = 40 mm
  else
    biasfwhm = opts.spm.biasfwhm;
  end
  if isempty(opts.spm.ngaus)
    % use default classes (a skull-stripped case is irrelevant here ;-)
    ngaus = boney_segment_preprocessing('ngaus',Ptpm);
  else
    ngaus = opts.spm.ngaus;
  end
  matlabbatch{1}.spm.spatial.preproc.channel.vols        = P;
  matlabbatch{1}.spm.spatial.preproc.channel.biasreg     = opts.spm.biasreg;
  matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm    = biasfwhm;
  matlabbatch{1}.spm.spatial.preproc.channel.write       = [0 1]; % write bias corrected image
  for ti = 1:numel(ngaus)
    matlabbatch{1}.spm.spatial.preproc.tissue(ti).tpm    = {sprintf('%s,%d',Ptpm,ti)};
    matlabbatch{1}.spm.spatial.preproc.tissue(ti).ngaus  = ngaus(ti);
    matlabbatch{1}.spm.spatial.preproc.tissue(ti).native = [ti<numel(ngaus) 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(ti).warped = [0 0];
  end
  matlabbatch{1}.spm.spatial.preproc.warp.mrf            = opts.spm.mrf;
  matlabbatch{1}.spm.spatial.preproc.warp.cleanup        = opts.spm.cleanup;
  matlabbatch{1}.spm.spatial.preproc.warp.reg            = opts.spm.warpreg;
  % subj was more robust compair to the default (=mri)
  matlabbatch{1}.spm.spatial.preproc.warp.affreg         = opts.spm.affreg;
  matlabbatch{1}.spm.spatial.preproc.warp.fwhm           = opts.spm.fwhm;
  % samp=5 was more robust in some UKB cases but more critical in other
  % datasets (e.g. ADNI), so we use the SPM default
  matlabbatch{1}.spm.spatial.preproc.warp.samp           = opts.spm.samp;
  matlabbatch{1}.spm.spatial.preproc.warp.write          = [0 0];
  matlabbatch{1}.spm.spatial.preproc.warp.vox            = NaN;
  matlabbatch{1}.spm.spatial.preproc.warp.bb             = [NaN NaN NaN; NaN NaN NaN];
end
%=======================================================================
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
  else
    matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr                    = 0.5 + 0.5*bias;
    matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr                     = 0.5;
  end
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.restypes.optimal = [1 0.3];
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.setCOM       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP          = 1070;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.affmod       = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr        = .5; % faster? - no, better robust 
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
      {fullfile(spm_file(which('cat12'),'fpath'), 'templates_MNI152NLin2009cAsym/T1.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.brainmask    = ...
      {fullfile(spm_file(which('cat12'),'fpath'), 'templates_MNI152NLin2009cAsym/brainmask.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.cat12atlas   = ...
      {fullfile(spm_file(which('cat12'),'fpath'), 'templates_MNI152NLin2009cAsym/cat.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.darteltpm    = ...
      {fullfile(spm_file(which('cat12'),'fpath'), 'templates_MNI152NLin2009cAsym/Template_1_Dartel.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shootingtpm  = ...
      {fullfile(spm_file(which('cat12'),'fpath'), 'templates_MNI152NLin2009cAsym/Template_0_GS.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regstr       = 0.5;
  else
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regmethod.shooting.shootingtpm = ...
      {fullfile(spm_file(which('cat12'),'fpath'), 'templates_MNI152NLin2009cAsym/Template_0_GS.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regmethod.shooting.regstr = 0.5;
  end
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.bb           = 12;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.vox          = 1.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtres            = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtmethod         = 'pbt2x';
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.SRP               = 30;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.vdist             = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.scale_cortex      = 0.7;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.add_parahipp      = 0.1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.close_parahipp    = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.experimental        = 0; 
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.new_release         = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy                = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors        = 1; %
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb                = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print               = 2;
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
