function [Pout,out] = boney_segment(job)
%boney_segment. Extract bone measures based on SPM segmentation.
% ... DETAILED DESCRIPTION ...
%
%  [Pout,out] = boney_segment(job);
%
%  job         .. SPM job structure
%   .files     .. raw images (or SPM m*-files)
%   .opts      .. structure of input values
%    .ctpm     .. TPM: 0-automatic (default), 1-SPM TPM, 2-children TPM,
%                 3-BlaiottaTPM, 4-optimized TPM7blr (if available)
%    .pmethod  .. preprocessing method: 1-SPM12, 2-CAT12 [, 3-CTseg]
%    .bmethod  .. method used for evaluation (default=2):
%                  0 - SPM seg8t mat values evaluation (fast)
%                  1 - volume-based 
%                  2 - volume- and surface-based 
%    .refine   .. refine SPM segmentation (default=1)
%    .prerun   .. use existing SPM/CAT preprocessing (0-ask if the SPM 
%                 parameters differ, 1-always rerun, 2-use existing, 
%                 3-rerun if the SPM parameters differ)
%    .rerun    .. rerun boney processing (default=0)
%    .nproc    .. number of parallel processes (default=0)
%    .verb     .. display progress
%                 (0 - be silent, 1 - one line per subject, 2 - details)
%    .affreg   .. do affine registration based (default=0)
%    .reduce   .. voxel binning for surface creation (higher values create
%                 surfaces with less vertices, i.e., details)
%                 - reduction factor vs. vertices in humans:
%                     1~120k, 2~30k, 3~13k, 4~7k, 6~3k, 8~2k
%                 - robust results for 1-4
%    .Patlas   .. atlas for regional measures 
%    .Pmask    .. mask problematic regions 
%    .spm      .. SPM preprocessing parameters (tpm, ngaus, biasreg, 
%                 biasfwhm, samp, mrf, cleanup, warpreg, affreg, fwhm)
%   .output    .. structure of output values
%    .report   .. create JPG report file (default=3)
%    .writevol .. write volume output (default=0)
%    .writesurf.. write surface data (default=0)
%    .writeseg .. keep SPM segmentation (default=0)
%
%  Pout         .. output files for SPM dependencies 
%  out          .. output structure
%   .spm8       .. SPM data from the Unified Segmentation
%   .tis        .. global estimated tissue values based on SPM
%   .tismri     .. MRI-based tissue values 
%   .vROI       .. volume-based global/regional evaluation of the bone/head
%   .sROI       .. surface-based global/regional evaluation of the bone/head
%   .main       .. most relevant values 
%   .qc         .. quick QC flags (see boney_segment_QC)
%   .tpmsel     .. automatic TPM selection (see boney_segment_selectTPM)
%
%
%  Structure by subfunctions [additional] marked by == NAME ==
%   (1)    boney_segment_filenames
%   (2)    boney_segment_checkPP
%   (3)    boney_segment_preprocessing [boney_segment_selectTPM]
%   (4)    boney_segment_prepare_print
%   (5)    boney_segment_get_segmat
%   (6.1)  boney_segment_loadMRI
%   (6.2) [boney_segment_simpleBone] (classic/prototype)
%   (6.3)  boney_segment_evalSPMseg
%   (6.4) [boney_segment_refineSPM]
%   (6.5)  boney_segment_extractbone
%   (7)   [boney_segment_create_bone_surface]
%   (8)    boney_segment_QC
%   (9)    boney_segment_cleanup
%   (10)   boney_segment_cmdline
%
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________
% $Id$

  global cat_err_res boned; %#ok<GVMIS> % for CAT error report

  if ~exist('job','var'), help boney_segment; return; end

  P = job.files;
  if isempty(P) || isempty(P{1})
    help boney_segment;
    return;
  end


  % use the global boney settings of the main process in parallel processes
  if isfield(job,'process_index') && isfield(job,'boned')
    boned = job.boned;
  end


  % RD202403: currently force expert setting for rerun/lazy settings
  oldexpert = cat_get_defaults('extopts.expertgui'); 
  cat_get_defaults('extopts.expertgui',2)


  % == default parameters ==
  %  (the defaults should be equal to those of the GUI in tbx_cfg_boney)
  def.files             = {};       % SPM m-files
  def.opts.verb         = 1;        % display progress (0 - be silent, 1 - one line per subject, 2 - details)
  def.opts.pmethod      = 1;        % preprocessing method: 1-SPM12, 2-CAT12 [, 3-CTseg]
  def.opts.bmethod      = 2;        % method: 0 - SPM seg8t eval only, 1 - volume-based, 2 - surface based
  def.opts.ctpm         = 0;        % TPM selector: 0 - automatic selection, 1 - default SPM TPM for adults,
                                    %               2 - children TPM, 3 - BlaiottaTPM, 4 - TPM7blr (if available)
  def.opts.nproc        = 0;        % number of parallel processes (0 - no parallel processing)
  def.opts.prerun       = 0;        % SPM/CAT preprocessing if previous results exist:
                                    %  0 - use them but ask (GUI) or warn (batch) if the SPM parameters differ,
                                    %  1 - always rerun preprocessing,
                                    %  2 - use them (warn if the SPM parameters differ),
                                    %  3 - rerun preprocessing if the SPM parameters differ
  def.opts.rerun        = 0;        % rerun processing (0 - no load previous result if available, 1 - yes)
  def.opts.affreg       = 0;        % do affine registration based on new skull-stripping
  def.opts.bias         = 1;        % strong bias correction
  def.opts.reduce       = 2;        % voxel binning for surface creation
                                    % (higher values create surfaces with less vertices, i.e., details)
                                    %  - reducion factor vs. vertices in humans:
                                    %    1~120k, 2~30k, 3~13k, 4~7k, 6~3k, 8~2k
                                    %  - robust results for 1-4
  def.opts.refine       = 1;        % refine segmenation (0 - no, 1 - yes)
  def.opts.bnorm        = 'WM';     % WM, CSF, muscle, bone, fat ..
  def.opts.normCT       = 0;        % use hard defined CT tissue thresholds for CTseg (used flag?) - RD20250116 - true is not working!
  def.opts.Patlas       = {fullfile(spm('dir'),'toolbox','boney','boney_KADA_bone-regions-ordered.nii')};
  def.opts.Pmask        = {fullfile(spm('dir'),'toolbox','boney','boney_KADA_bone-mask.nii')};
  def.opts.snspace      = [80,7,2];                      % cmd print format - only internal
  def.opts.MarkColor    = cat_io_colormaps('marks+',40); % colormap for ratings - only internal
  def.opts.reslim       = 1.5;      % general resolution limit
  def.opts.expert       = 2;        % user level (0 - default, 1 - expert, 2 - developer)
  def.opts.classic      = 0;        % estimate also first prototype version
  % SPM preprocessing parameters (empty values are defined by the TPM/opts.bias)
  def.opts.spm.tpm      = '';       % TPM file (overwrites opts.ctpm)
  def.opts.spm.ngaus    = [];       % number of Gaussians per class (default [1 1 2 3 4 2] and [1 1 1 3 1 2 3] for the BlaiottaTPM)
  def.opts.spm.biasreg  = 0.001;    % bias regularisation
  def.opts.spm.biasfwhm = [];       % bias FWHM (default: 60 - 20*opts.bias)
  def.opts.spm.samp     = 3;        % sampling distance
  def.opts.spm.mrf      = 1;        % MRF cleanup
  def.opts.spm.cleanup  = 1;        % light cleanup
  def.opts.spm.warpreg  = [0 0.001 0.5 0.05 0.2]; % warping regularisation
  def.opts.spm.affreg   = 'subj';   % affine regularisation
  def.opts.spm.fwhm     = 0;        % smoothness
  def.output.report     = 3;        % write report
  def.output.resdir     = 'BIDS';   % result directory
  def.output.writevol   = 0;        % write volume output
  def.output.writeseg   = 0;        % keep SPM/CAT segmentation (0 - no, 1 - yes, 2 - as mat file)
  def.output.writesurf  = 0;        % write surface data
  def.output.writeCSV   = 1;
  def.output.prefix     = 'boney_';
  job                   = cat_io_checkinopt(job,def);
  if ~isfield(job,'ppfallback'), job.ppfallback = false; end % internal: segmentation of a TPM fallback
  job.opts.spm          = cat_io_checkinopt(job.opts.spm,def.opts.spm);
  job.opts.subdirs      = job.opts.pmethod == 2 && cat_get_defaults('extopts.subfolders');
  job.output.resdir     = strrep(job.output.resdir,'BIDS',['..' filesep 'derivatives' filesep 'boney']);
  job.output.resdir     = strrep(job.output.resdir,'/',filesep);

  Pout = struct();

  % no surface output (report==3) without surface processing (bmethod == 2)
  job.output.report     = min(job.output.report,max(2,job.opts.bmethod + 1));

  % filenames for dependencies
  % job.opts.fmethod .. method defined by the selected file
  % (files that cannot be processed are removed)
  [out,job.opts.fmethod,job.opts.pmethod,P] = boney_segment_filenames(P,job);
  job.files = P;

  % == check existing SPM preprocessing ==
  %  Test if SPM/CAT preprocessing is required or if existing segmentations
  %  can be used. Differences between the SPM parameters of the existing
  %  segmentation and the current setting are reported and the user can
  %  decide how to proceed (only in the main process).
  [pp,job,stopit] = boney_segment_checkPP(out,job);
  if stopit
    cat_get_defaults('extopts.expertgui',oldexpert)
    return
  end
  pproc = pp.proc; 

  % cases that will be processed (and not just loaded) in the main loop below
  % i.e. no XML, XML older than the image or the segmentation, or new segmentation
  rproc = true(1,numel(P));
  if ~job.opts.rerun
    for i = 1:numel(P)
      rproc(i) = pp.force(i) || boney_isolder(out(i).P.xml, {out(i).P.org, out(i).P.seg8}); 
    end
  end


  % == parallel processing ==
  %  Split the files that have to be processed into separate MATLAB processes
  %  (mostly to run the SPM segmentation in parallel) and reload all results
  %  afterwards in this process to write one CSV file and to return the
  %  output for dependencies.
  if job.opts.nproc > 0 && ~isfield(job,'process_index') && sum(rproc) > 1
    pjob       = job;
    pjob.files = P(rproc);
    pjob.nproc = job.opts.nproc;
    pjob.boned = boned;  % the global boney settings are not transferred by cat_parallelize
    cat_parallelize(pjob,mfilename,'files');

    % reload all results (or process failed cases) but use the existing
    % segmentations (just created by the parallel processes)
    job.opts.nproc  = 0;
    job.opts.rerun  = 0;
    job.opts.prerun = 2;
    cat_get_defaults('extopts.expertgui',oldexpert)
    [Pout,out] = boney_segment(job);
    return
  end


  if job.opts.fmethod && any((pproc | pp.bconly) & rproc)
    % call SPM/CAT preprocessing if raw images are used but only for cases
    % that will be processed (and not just loaded) in the main loop below
    % (if only the bias corrected image is missing, it is written by the 
    % existing SPM parameters)
    boney_segment_preprocessing(P(pproc & rproc), out(pproc & rproc), job, out(pp.bconly & rproc));

    % update filenames (e.g. the number of classes depends on the TPM)
    out = boney_segment_filenames(P,job);
  end

  % =======================================================================
  % TODO:
  % * add real fat measures (first version of head fat fraction and fat 
  %   thickness in boney_segment_extractbone needs validation)
  % * use intnormed image? (intra- vs. intersite variability)
  % * use fat segmentation to avoid critical regions
  % * further postprocessing in case of children (use children TPM as
  %   indicator) to correct the head-scull misalignment
  % * image boundaries (e.g. IBSR 18), defacing, and skull-stripping are 
  %   only detected by boney_segment_QC but not corrected
  % * log/exp intensity normalization in "boney_segment_loadMRI" 
  % * report: colorbar, vROI/sROI table, TIV, line-plots
  % * CT (CTseg) testing and debugging 
  % * normalize bone/head measures, bone-relation measure, fat-relation measure
  % * estimation of fat-shift to quantify understimation of bone ?
  % * maybe use the SPM bone marrow segmentation error as feature to 
  %   classify the fat tissue more directly, e.g. with 
  %     Gauss = [  1  1  1-2  1-2 3-4  2  ]
  %   but this would not allow to use the easy SPM defaults
  % =======================================================================

  %#ok<*CLOCK,*DETIM>

  % == create table elements for the in-line report ==
  %  - just give some short overview about the processing request and parameters
  %  - Feature (visual):  one line progress >> first name, then progress + processing dots/stars ...
  %  - Feature (expert):  print additional parameters
  Theader = boney_segment_prepare_print(P,job);
  stime3 = clock; i = 1; %#ok<*NASGU>



  %% main loop over all subjects
  spm_progress_bar('Init',numel(job.files),'Bone Segment','Volumes Complete');
  for i = 1:numel(P)
    % reset global CAT error variable
    clearvars -global cat_err_res;
    stime2                   = clock;
    cat_err_res.stime        = clock;
    cat_err_res.cat_warnings = cat_io_addwarning('reset'); % reset warnings

    % remove subject specific results of the previous subject
    clear vROI sROI Si Stm Ye Vo Yo Ym Yc Ya Ymsk Ybonemarrow tismri Affine seg8t tis;


    if rproc(i)
     try
      % == GET SPM DATA ==
      %  - get and evaluate the original SPM preprocessing structure (seg8t)
      %    extract further values (tis) and voxel size
      stime = cat_io_cmd('  Load SPM','g5','',job.opts.verb>1);
      [ seg8t , tis , vx_vol , trans ] = boney_segment_get_segmat( out(i), job.output.writevol, job.opts.normCT, job.opts.verb );
      if isempty(vx_vol) || numel(fieldnames(tis))==0
        error('boney:segment:noSPMdata','Cannot use the SPM/CAT preprocessing data.');
      end

      % quick QC of the SPM preprocessing based on the seg8t values
      tis = boney_segment_QC('seg8',tis,seg8t);

      % subject specific settings
      jobi = job;
      if out(i).CTseg || tis.weighting<0
        jobi.opts.affreg = -1; % CT data require a new affine registration
      end
      if seg8t.ntpm == 7 % BlaiottaTPM refinement requires the SPM TPM bone class as second mask
        jobi.opts.Pmask{2} = fullfile(spm('dir'),'tpm','TPM.nii,4');
      end


      % == major branch for bone method ==
      %  0 only SPM values (i.e., super fast)
      %  1 use MRI
      %  2 use MRI + create surfaces (i.e., a bit slower but much cooler although it does not support SBM yet)
      if job.opts.bmethod


        % == load MRI and segmentation ==
        %  - load the (resolution-limited) images
        %  - show images for debugging:
        %      ds('d2sm','',vx_vol,Ym .* (0.5 + 0.5*Ybone) ,Ya,90)
        stime = cat_io_cmd('  Load MRIs','g5','',job.opts.verb>1,stime);
        [Vo,Yo,Yc,Ya,Ymsk,Ym, Affine, YaROIname, YaROIid, RES, BB] = ...
          boney_segment_loadMRI( out(i).P, jobi, seg8t, tis, 1:numel(out(i).P.cls), 25); % CT rand issue
        spm_progress_bar('Set',i - 0.8); vx_vol = RES.vx_volr;


        % == BlaiottaTPM ==
        %  - reduce the 7 BlaiottaTPM classes to the class 6 TPM classes
        if seg8t.ntpm == 7
          bb=2; sY = size(Yc{1}); Ybb = true(size(Yc{1})); Ybb(bb+1:sY(1)-bb,bb+1:sY(2)-bb,bb+1:sY(3)-bb) = false; clear bb sY; 
          Ybg   = smooth3( max( Yc{7} , cat_vol_morph(cat_vol_morph(Yc{7}>.2 | (Ym<.5 & Ybb),'dc',2),'de',2) )).^2; clear Ybb; 
          Yc{4} = Yc{6} .* max( cat_vol_morph(Ymsk{2}>.1,'dd',2) & Ym<.5, Ymsk{2}>.1)  + ... hard bone
                  (Yc{5}+Yc{4}) .* smooth3(cat_vol_morph(Ymsk{2}>.5 & Ym>.3).^.5); % bone marrow
          Yc{6} = single(Ybg);
          Yc{5} = max(0, 1 - (Yc{1}+Yc{2}+Yc{3}+Yc{4}+Yc{6}) );
          Yc(7) = [];
          % the Gaussians (seg8t.lkp) were already mapped in boney_segment_get_segmat
        end

        % == original measure from our prototype concept ==
        if job.opts.classic
          stime = cat_io_cmd('  Classic bone measures','g5','',job.opts.verb>1,stime);
          out(i).classic = boney_segment_simpleBone(seg8t,Yo,Yc,job.opts.refine,vx_vol); 
          spm_progress_bar('Set',i - 0.75);
        end


        % == evaluate/test the SPM segmentation ==
        %  - full estimation with refined tissue peaks
        stime  = cat_io_cmd('  Evaluate MRIs','g5','',job.opts.verb>1,stime);
        [tismri, Ybraindist0] = boney_segment_evalSPMseg(Yo,Ym,Yc,Ymsk{1},vx_vol,0,job,seg8t,tis);
        spm_progress_bar('Set',i - 0.7);


        % == refine SPM segmentation or just prepare some maps ==
        %  - store the changes also in the output structure
        if job.opts.refine && ~out(i).CTseg && tis.weighting >= 0 % ~CT
          if all(tis.res_vx_vol < 1.5)
            stime = cat_io_cmd('  Refine SPM','g5','',job.opts.verb>1,stime);
            [Yc,Ye,Ya,out(i).boney_refine_SPM] = boney_segment_refineSPM(Yo,Ym,Yc,Ya,Ybraindist0,tis,tismri,vx_vol);
          else
            stime = cat_io_cmd('  No refinement because of low (slice) resolution!','g5','',job.opts.verb>1,stime);
            cat_io_addwarning('Warning:TooLowResForRef','No refinement because of too low resolutions!',1,[1 1],{},0,job.opts.verb>1);
            Ye = cell(0);
          end
        else
          Ye = cell(0);
        end
        spm_progress_bar('Set',i - 0.6);
        clear Yo

        
        % == get bone measures ==
        %  - get different bone maps Y* and the volume-based regional values vROI:
        %      Y*   .. bone/head maps for surface mapping
        %      vROI .. extracted global/regional bone/head values
        stime = cat_io_cmd('  Extract bone measures','g5','',job.opts.verb>1,stime);
        [Ybonepp,Ybonethick,Ybonemarrow,Yheadthick, vROI] = ...
          boney_segment_extractbone(Vo,Ym,Yc,Ye,Ya,Ymsk{1},trans,seg8t,tis,tismri,out(i),jobi,vx_vol,YaROIname,YaROIid,RES,BB);
        spm_progress_bar('Set',i - 0.4);


        % == restore resolution & boundary box ==
        Ym            = cat_vol_resize( Ym            ,'dereduceV'    ,RES);
        Ymsk{1}       = cat_vol_resize( Ymsk{1}       ,'dereduceV'    ,RES,'nearest');
        [Ym, Ymsk{1}] = cat_vol_resize( {Ym, Ymsk{1}} ,'dereduceBrain',BB);
        for ai = 1:numel(Ya)
          Ya{ai} = cat_vol_resize( Ya{ai}  ,'dereduceV'    ,RES,'nearest'); % labels
          Ya{ai} = cat_vol_resize( Ya{ai}  ,'dereduceBrain',BB);
        end
        for ci = 1:numel(Yc)
          Yc{ci} = cat_vol_resize( Yc{ci}  ,'dereduceV'    ,RES);
          Yc{ci} = cat_vol_resize( Yc{ci}  ,'dereduceBrain',BB);
        end


        % == surface-based processing ==
        %  - this function creates the bone surface with different measures
        %    (two surface S* are used for the report) and the surface-based
        %    regional values sROI:
        %      Si/Stm .. bone-surface with bone-intensity/bone-thickness
        %      sROI   .. extracted global/regional bone/head surface values
        %  - although the surfaces can be saved, it does not suport surface
        %    registration now!
        %  - only useful for not too low resolutions, i.e. slice res < 2.1 and slice thickness < 5.1
        if job.opts.bmethod>1  
          if  sum(vx_vol <= 2.1)>1  &&  all(vx_vol < 5.1) 
            stime = cat_io_cmd('  Extract bone surfaces','g5','',job.opts.verb>1,stime);
            [Si, Stm, sROI] = boney_segment_create_bone_surface(Vo, Ybonepp, Ybonemarrow, Ybonethick, Yheadthick, Ya, Ymsk{1}, YaROIname, YaROIid, out(i), jobi);
          else
            Si = ''; Stm = '';
            cat_io_addwarning('Warning:TooLowResForSurf', ['No surface processing because of too low resolutions, i.e. \\n' ...
               'slice resolution < 2.1 and slice thickness < 5.1, use volume measures!'],2,[1 1],{},0,job.opts.verb>1);
          end
        else
          Si = ''; Stm = '';
        end

      else
        % == fast SPM mat-based pipeline ==
        %  - fast version that only looks for (and exports) the bone tissue
        %    values and creates some report (Yo==Ym?)
        %  - no other tissues (e.g. GM etc.) are loaded and processed!
        [Vo, Ym, Yc, Ybonemarrow, tismri, Si, Stm, Affine] = ...
          boney_segment_fst(out(i).P, jobi, seg8t, tis, vx_vol);
      end
      spm_progress_bar('Set',i - 0.2);



      % == create output structure
      if isfield(seg8t.image,'private'), seg8t.image = rmfield(seg8t.image,'private'); end
      if isfield(seg8t,'tpm') % SPM but not CAT
        if isfield(seg8t.tpm,  'private'), seg8t.tpm   = rmfield(seg8t.tpm  ,'private'); end
      end
      if isfield(seg8t,'tpmA'), seg8t = rmfield(seg8t,'tpmA'); end
      if isfield(seg8t,'sett'), seg8t = rmfield(seg8t,'sett'); end
      out(i).spm8    = seg8t;
      out(i).tis     = tis;
      out(i).tismri  = tismri;
      out(i).opts    = job.opts;
      out(i).date    = datestr(stime3,'YYYYmmDD-HHMMSS'); %#ok<*DATST>
      if exist(out(i).P.tpmsel,'file') % automatic TPM selection 
        out(i).tpmsel = load(out(i).P.tpmsel); 
      end
      if exist('vROI','var'), out(i).vROI = vROI;   end
      if exist('sROI','var'), out(i).sROI = sROI;   end


      % == final measures for the bone and head ==
      %  * [v|s]BMDH have to use the inverted measure as the increased  
      %    intensity of fatty marrow results in weaker bones (lower BMD) 
      %    with H for the UKB BMD head measure
      if exist('vROI','var') && numel(out(i).vROI(1).bonecortex)>2
        out(i).main.vBMDH = -out(i).vROI(1).bonecortex(3);
      elseif exist('vROI','var')
        out(i).main.vBMDH = -out(i).vROI(1).bonecortex(end);
      end
      if exist('sROI','var') && numel(out(i).sROI(1).bonecortex)>2
        out(i).main.sBMDH = -out(i).sROI(1).bonecortex(3);
      elseif exist('sROI','var')
        out(i).main.sBMDH = -out(i).sROI(1).bonecortex(end);
      end

      % == quality control ==
      %  quick QC flags of the SPM preprocessing and the bone/head measures
      %  (used for warnings and the selection of the TPM)
      out(i).qc = boney_segment_QC('final',out(i)); 
      out(i).ptime = etime(clock,stime2); % processing time (without SPM/CAT preprocessing)


      % == automatic TPM fallback ==
      %  A failed segmentation (see boney_segment_QC) of the automatic TPM 
      %  selection is processed again with the next TPM (BlaiottaTPM, 
      %  TPM7blr if available, SPM TPM) until the QC passes or all TPMs were
      %  tested (the last result is kept). Only segmentations created by 
      %  boney are replaced. The XML and the report are written by the 
      %  recursive call. 
      Ptpmfb = '';
      if isfield(out(i).qc,'segfail') && out(i).qc.segfail && isempty(job.opts.spm.tpm) && ...
         job.opts.ctpm == 0 && job.opts.pmethod == 1 && (pproc(i) || job.ppfallback) && ...
         isfield(out(i),'tpmsel') && isfield(out(i).tpmsel,'ctpm')
        [Ptpmfb, tpmsel] = nextTPM(out(i).tpmsel);
      end
      if ~isempty(Ptpmfb)
        if job.opts.verb
          cat_io_cprintf('warn',sprintf('  QC failed (%s) - rerun SPM with %s.\n', ...
            strjoin(out(i).qc.flags,', '), spm_str_manip(Ptpmfb,'t')));
        end
        boney_segment_preprocessing(P(i), out(i), job, [], Ptpmfb);
        save(out(i).P.tpmsel,'-struct','tpmsel');

        optb                 = job; 
        optb.files           = job.files(i);
        optb.opts.verb       = 0; 
        optb.opts.rerun      = 1; 
        optb.opts.nproc      = 0; 
        optb.opts.prerun     = 2;    % use the new segmentation
        optb.output.writeCSV = 0;
        optb.ppfallback      = true; % segmentation created by boney
        [~,outi]             = boney_segment(optb);
        out                  = cat_io_mergeStruct(out,outi,[],i);
        rerunstr             = sprintf('fallback %s',spm_str_manip(Ptpmfb,'tr'));
      else
        % export to XML
        %  An error of another subject adds an (empty) error field to all
        %  elements of the out struct array that has to be removed here.
        if isfield(out,'error'), out(i).error = []; end
        outi = out(i); if isfield(outi,'error'), outi = rmfield(outi,'error'); end
        cat_io_xml(out(i).P.xml, outi); % export to XML


        %% == create report ==
        if job.output.report
          stime = cat_io_cmd('  Create report','g5','',job.opts.verb>1,stime);
          boney_segment_print_figure(Vo, Ym, Yc, Ybonemarrow, Si, Stm, out(i), job, Affine);
        end
        if job.opts.verb>1, fprintf('% 5.0fs \n',etime(clock,stime)); end 
        spm_progress_bar('Set',i - 0.1);
        rerunstr = '';
      end

     catch e
      % == processing error ==
      %  Write the error into the XML to mark the case as processed but 
      %  failed (the case is reprocessed only by rerun=1).
      out(i).error = struct('message',e.message,'identifier',e.identifier, ...
        'stack',strjoin(arrayfun(@(s) sprintf('%s:%d',s.name,s.line),e.stack,'UniformOutput',false),' < '));
      out(i).date  = datestr(stime3,'YYYYmmDD-HHMMSS');
      cat_io_xml(out(i).P.xml, struct('P',out(i).P,'CTseg',out(i).CTseg,'error',out(i).error,'date',out(i).date)); 
      rerunstr = 'failed';
     end

    else
    % == load previous result, reprocess in case of error ==
      try
        outi     = cat_io_xml(out(i).P.xml);
        rerunstr = 'loaded';
        if isfield(outi,'error') && ~isempty(outi.error), rerunstr = 'failed'; end
      catch
        optb              = job; 
        optb.files        = job.files(i);
        optb.opts.verb    = 0; 
        optb.opts.rerun   = 1; 
        optb.opts.nproc   = 0; 
        optb.output.writeCSV = 0;
        [~,outi]          = boney_segment(optb);
        rerunstr          = 'rerun';
      end
      out    = cat_io_mergeStruct(out,outi,[],i);
    end


    % == prepare output and cleanup files ==
    % (the segmentation is only removed if it was created here) 
    Pout = boney_segment_cleanup(Pout,out,job,i,pproc(i) && rproc(i));


    % == print command line report ==
    if isfield(out(i),'error') && ~isempty(out(i).error)
      if job.opts.verb
        cat_io_cprintf('err',sprintf('%4d) %s: %s\n',i, ...
          spm_str_manip(job.files{i},['a' num2str(job.opts.snspace(1) - 14)]),out(i).error.message));
      end
    else
      boney_segment_cmdline(job,out,i,stime2,rerunstr);
    end
    spm_progress_bar('Set',i);
    
  end

  % final print and cleanup
  if job.opts.verb > 0
    fprintf('%s\n',repmat('-',size(Theader)));
    fprintf('Boney main processing done.\n\n'); 
  end
  spm_progress_bar('Clear');



  % == final CSV export of all subjects ==
  %  (not in the parallel processes, as the main process writes one CSV for 
  %  all, and only for successfully processed cases)
  if job.output.writeCSV && ~isfield(job,'process_index') && isfield(Pout,'xml') 
    Pcsv = {}; 
    for i = 1:min(numel(out),numel(Pout.xml))
      if ~isempty(Pout.xml{i}) && ~( isfield(out(i),'error') && ~isempty(out(i).error) )
        Pcsv{end+1} = Pout.xml{i}; %#ok<AGROW>
      end
    end
    if numel(Pcsv)>1
      if job.opts.verb
        stime = cat_io_cmd('Create CSV with all subjects','g9','',job.opts.verb>0);
      end
      rlevel = {'boney_default','boney_expert'};
      matlabbatch{1}.spm.tools.boney.xml2csv.files        = Pcsv';
      matlabbatch{1}.spm.tools.boney.xml2csv.outdir       = {out(1).P.resdir}; % result directory (not the raw data)
      matlabbatch{1}.spm.tools.boney.xml2csv.fname        = sprintf('Boney_xml2csv_report_%s.csv',datestr(stime3,'YYYYmmDD-HHMMSS'));
      matlabbatch{1}.spm.tools.boney.xml2csv.fieldnames   = {' '};
      matlabbatch{1}.spm.tools.boney.xml2csv.avoidfields  = {''};
      matlabbatch{1}.spm.tools.boney.xml2csv.report       = rlevel{ ( boned.expertgui > 0) + 1 };

      % run silently
      try
        evalc('spm_jobman(''run'',matlabbatch);');
      catch
        cat_io_cprintf('err','Failed to write final CSV file! Check boney_xml2csv.m function.\n')
      end
      if job.opts.verb > 0, fprintf('% 5.0fs\n\n',etime(clock,stime)); end
    end
  end

  cat_get_defaults('extopts.expertgui',oldexpert)

end
function [Ptpm,tpmsel] = nextTPM(tpmsel)
%nextTPM. Next TPM of the automatic fallback that was not tested before 
% (see boney_segment_selectTPM, the TPM7blr is only used if available).
  Ptpms = boney_segment_preprocessing('tpms');
  if ~isfield(tpmsel,'tried'),    tpmsel.tried    = tpmsel.ctpm; end
  if ~isfield(tpmsel,'fallback'), tpmsel.fallback = 0; end
  Ptpm  = '';
  for ti = [3 4 1]
    if ~any(tpmsel.tried == ti) && exist(Ptpms{ti},'file')
      Ptpm = Ptpms{ti};
      if tpmsel.fallback == 0, tpmsel.fallback = tpmsel.ctpm; end % originally selected TPM
      tpmsel.ctpm  = ti;
      tpmsel.tried = [tpmsel.tried ti];
      return
    end
  end
end
function isold = boney_isolder(file,reffiles)
%boney_isolder. True if the file does not exist or is older than one of the 
% existing reference files (independent of the CAT expert level).
  isold = ~exist(file,'file');
  if ~isold
    fd = dir(file); 
    for ri = 1:numel(reffiles)
      if exist(reffiles{ri},'file')
        rd    = dir(reffiles{ri}); 
        isold = isold || fd.datenum < rd.datenum; 
      end
    end
  end
end
