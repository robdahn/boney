function [Pout,out] = boney_segment(job)
%boney_segment. Extract bone measures based on SPM segmentation.
% ... DETAILED DESCRIPTION ...
% 
%  Pout = boney_segment(job); 
%
%  job         .. SPM job structure
%   .files     .. SPM m-files
%   .opts      .. structure of input values
%    .bmethod  .. [SPM,CAT]
%    .bmethod  .. method used for evaluation (default=2): 
%                  0 - SPM seg8t mat values evaluation (fast)
%                  1 - SPM volumes with problems for high intensity bone 
%                      marrow 
%                  2 - refined volumes  
%    .verb      .. display progress 
%                 (0 - be silent, 1 - one line per subject, 2 - details)
%    .affreg    .. do affine registration based (default=1) 
%    .reduce    .. voxel binning for surface creation (higher values create 
%                 surfaces with less vertices, i.e., details)
%                 - reduction factor vs. vertices in humans:   
%                     1~120k, 2~30k, 3~13k, 4~7k, 6~3k, 8~2k
%                 - robust results for 1-4
%    .mask      .. mask problematic regions (default=1)
%   .output     .. structure of output values
%    .report    .. create JPG report file (default=1) 
%    .writevol  .. write volume output (default=0) 
%    .writesurf .. write surface data (default=0)
%    
%  out          .. output structure
%   .raw        .. original results (e.g., min, max, median, std ...) for the 
%                  different regions etc. 
%   .spm8       .. SPM data from the Unified Segmentation  
%   .tis        .. global estimated 
%   .vol        .. volume-based global/regional evaluation of the bone/head
%   .surf       .. surface-based global/regional evaluation of the bone/head
%   .main       .. most relevant values (see ...)
% 
%
%  Structure by subfunctions [additional] marked by == NAME ==
%   (1)    boney_segment_filenames
%   (2)    boney_segment_preprocessing
%   (3)    boney_segment_prepare_print
%   (4)    boney_segment_get_segmat
%   (5.1)  boney_segment_loadMRI
%   (5.2) [boney_segment_simpleBone] (classic/prototype)
%   (5.3)  boney_segment_evalSPMseg
%   (5.4) [boney_segment_refineSPM]
%   (5.5)  boney_segment_extractbone
%   (6)   [boney_segment_create_bone_surface]
%   (7)    boney_segment_cleanup
%   (8)    boney_segment_cmdline
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


  % == default parameters ==
  def.files             = {};       % SPM m-files
  def.opts.verb         = 2;        % display progress (0 - be silent, 1 - one line per subject, 2 - details)
  def.opts.pmethod      = 1;        % preprocessing method: 1-SPM12, 2-CAT12 [, 3-CTseg] 
  def.opts.bmethod      = 2;        % method: 0 - SPM seg8t eval only, 1 - volume-based, 2 - surface based
  def.opts.ctpm         = 1;        % TPM selector: 1 - default TPM for adults; 2 - children TPM; 
  def.opts.prerun       = 0;        % avoid preprocessing in matching files of SPM/CAT are available
  def.opts.rerun        = 0;        % rerun processing (0 - no load previous result if available, 1 - yes) 
  def.opts.affreg       = 0;        % do affine registration based on new skull-stripping  
  def.opts.bias         = 1;        % strong bias correction
  def.opts.reduce       = 4;        % voxel binning for surface creation 
                                    % (higher values create surfaces with less vertices, i.e., details)
                                    %  - reducion factor vs. vertices in humans:   
                                    %    1~120k, 2~30k, 3~13k, 4~7k, 6~3k, 8~2k
                                    %  - robust results for 1-4
  def.opts.refine       = 1;        % refine segmenation 
  def.opts.subdirs      = 1;        % use subdirectories 
  def.opts.normCT       = 0;        % use hard defined CT tissue thresholds for CTseg (used flag?) 
  def.opts.Patlas       = {fullfile(spm('dir'),'toolbox','boney','boney_KADA_bone-regions.nii')};
  def.opts.Pmask        = {fullfile(spm('dir'),'toolbox','boney','boney_KADA_bone-mask.nii')};
  def.opts.snspace      = [80,7,2];                      % cmd print format - only internal
  def.opts.MarkColor    = cat_io_colormaps('marks+',40); % colormap for ratings - only internal
  def.opts.reslim       = 1.5;      % general resolution limit 
  def.opts.expert       = 2;        % user level (0 - default, 1 - expert, 2 - developer)
  def.opts.classic      = 1;        % estimate also first prototype version 
  def.output.report     = 2;        % write report
  def.output.writevol   = 1;        % write volume output 
  def.output.writeseg   = 1;        % write volume output 
  def.output.writesurf  = 1;        % write surface data 
  def.output.writeCSV   = 1; 
  def.output.prefix     = 'boney_';
  job                   = cat_io_checkinopt(job,def);
  Pout = struct();

  % no surface output (report==3) without surface processing (bmethod == 2)
  job.output.report     = min(job.output.report,max(2,job.opts.bmethod + 1)); 

  % filenames for dependencies 
  [out,job.opts.fmethod,job.opts.pmethod] = boney_segment_filenames(P,job);


% ######################  
% just some fast test settings
if 0
  job.opts.verb    = 1;
  job.opts.classic = 1;
  job.opts.rerun   = 0; 
end
% ######################  



  if job.opts.fmethod 
    % call SPM/CAT preprocessing if raw images are used  
    boney_segment_preprocessing(P, out, job.opts.ctpm, job.opts.pmethod, job.opts.bias, job.opts.prerun);
  end



   
  % =======================================================================
  % TODO: 
  % * Check volumes and CSV to detect bad images (QC?) and segmentations 
  %   (soft values and high variance in values)
  %   >> warning + report + suggestion (affine, resolution)
  %
  % * add real fat measures 
  %   - bone vs. head thickness? no, you will need fat vs. skin/muscles
  %   - required for better scaling at some point ...
  %
  % * use intnormed image?
  %   - our intensity values are going over different classes - intra- vs. intersite variability 
  %   - use CSF-WM contrast (but what about HH?) 
  %
  % * use fat segmentation to avoid critical regions
  %
  % * test SPM progress bar
  % * add rerun parameter
  %
  % * test Brudfors segmentation 1 and 2  
  % * update Yc map (at least print this one) ... diff muscles/bone?
  % * Ym-scaling in MT incrorrect
  % * colorbar vol from reportfig !!!!!!!!!!
  % * colorbar histogram for surf
  % * simple version just to report SPM segmentation (=simplest case)
  %   > render skull + brain with pbt 1 mm with thicknes ;-)
  %
  % * list of used CAT functions ?
  % 
  % * further postprocessing in case of children (use children TPM as 
  %   indicator) to correct the head-scull misalignment
  %
  %  * bone / bone marrow TPM concept 
  %    > this was not really working and it is probably better to move on
  %  * children case concept
  %    > done by cat_children TPM but it would be nice to have some younger TPM 
  %  * handle image boundaries by NaNs eg. IBSR 18
  %    > this works more or less
  %    > a warning would be good 
  %      - detection based on the BB? and Affine registration in "boney_segment_get_segmat" ?
  %      - detection in volumes i.e. in "boney_segment_loadMRI" or following functions ?
  %  * handle defaced areas 
  %    > detect defacing "boney_segment_loadMRI"
  %  * catch skull-stripping case (skull-stripping / deface defintions) 
  %    > if the head class contains too much background like values the TPM 
  %      is not fitting > could be used for detection and maybe correction
  %      (but probably more detection and suggestion of other TPM)
  %      - should be tried in "boney_segment_get_segmat"
  %    > detect defacing "boney_segment_loadMRI"
  %
  %  * write table as csv and add conclusion on command line
  %    > call xml2csv batch internally? If so than which level? > GUI level
  %
  %  * log/exp intensity normalization in "boney_segment_loadMRI" ##########
  %    
  %  * print vROI/sROI table in report
  %  * add TIV in figure?
  %  * add line-plots for special bone/head measures ? or maybe as extra signs + 
  %    maybe only print additional peaks on top ? 
  %    normalize each tissue peak by its volume ?
  %
  %  * Public (one MRI and one CT from OASIS with 1.5 mm) + 
  %    private test data (UKB subsample (long) + CT-OASIS-long)  
  %  * Test-scripts for parameters (just run) and documentary of major 
  %    parameters with public/private results/analysis.
  %  * estimation of fat-shift to quantify understimation of bone ?
  %  * add BIDS !
  %
  % ##############################################
  %  * CT Testen + Debuggen ################
  %  * figure output vROI and sROI+ classic #############
  %  * cmd output vROI and sROI ################
  %  * normalize SPM bone/head measures
  %  * normalize bonemed / sROI / vROI measures??? 
  %  * bone-relation measure, fat-relation measure, relative volumes ... thickness 
  % 
  % =======================================================================

  % =======================================================================
  % ##### cleanup #######
  % FEATURES are my BUG :(
  %
  % What do I realy need!
  % - fast mat-based intensity measurements [ SPM8MAT > TIS ]
  %    - global kmeans
  % - (refined) volumen-based intensity measurements [ TISMRI ] 
  %    - global simple (median)
  %    - global refined (sub-classes by kmeans) 
  %    - global/regional simple mean
  %    - histogram for ML?
  % - volume and thickness measurements of bone[, head, and brain] 
  %
  % CASES: ###########################################
  % - superfast
  %   - only mat values 
  %   - no img Report
  % - fast global
  %   - fast bone (median intensity, lowres thickness in limited range, kmeans2)
  %   - fast head (median intensity, lowres thickness in limited range, kmeans3)
  %   - partial img report
  % - refined segments + regions
  %   - refined bone (median intensity?, lowres thickness in limited range, kmeans2)
  %   - refined head (median intensity?, lowres thickness in limited range, kmeans3)
  %   - full img report
  % - surface based (avg,min,max)
  %   
  % =======================================================================


  % == create table elements for the in-line report == 
  %  - just give some short overview about the processing request and parameters 
  %  - Feature (visual):  one line progress >> first name, then progress + processing dots/stars ... 
  %  - Feature (expert):  print additional parameters
  boney_segment_prepare_print(P,job);
  stime3 = clock; i = 1; %#ok<NASGU> 
  


  %% main loop over all subjects
  spm_progress_bar('Init',numel(job.files),'Bone Segment','Volumes Complete');
  for i = 1:numel(P)
    % reset global CAT error variable 
    clearvars -global cat_err_res;
    stime2                   = clock;   
    cat_err_res.stime        = clock; 
    cat_err_res.cat_warnings = cat_io_addwarning('reset'); % reset warnings 


    if ... %cat_io_rerun(which(mfilename),out(i).P.xml) || ...
        cat_io_rerun(out(i).P.org,out(i).P.xml,0,-1) || job.opts.rerun
  
  
      % == GET SPM DATA ==
      %  - get and evaluate the original SPM preprocessing structure (seg8t)
      %    extract further values (tis) and voxel size 
      stime = cat_io_cmd('  Load SPM','g5','',job.opts.verb>1); 
      [ seg8t , tis , vx_vol ] = boney_segment_get_segmat( out(i), job.opts.verb ); 
      job.reslim = max(job.opts.reslim,mean(vx_vol)); % limit processing resolution 

  
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
if out(i).CTseg, job.affreg = -1; end % this is not optimal here - replace it later 
        [Vo,Yo,Yc,Ya,Ymsk,Ym, Affine, YaROIname, RES, BB] = ...
          boney_segment_loadMRI( out(i).P, job, seg8t, tis, 1:5, 25); % CT rand issue
        spm_progress_bar('Set',i - 0.8); vx_vol = RES.vx_volr; 
      

        % == original measure from our prototype concept ==
        if job.opts.classic
          stime = cat_io_cmd('  Classic bone measures','g5','',job.opts.verb>1,stime); 
          out(i).classic = boney_segment_simpleBone(seg8t,Yo,Yc); % ##### use Ym ?  #########
          spm_progress_bar('Set',i - 0.75);
        end

  
        % == evaluate/test the SPM segmentation == 
        %  - full estimation with refined tissue peaks
        stime  = cat_io_cmd('  Evaluate MRIs','g5','',job.opts.verb>1,stime); 
        [tismri, Ybraindist0] = boney_segment_evalSPMseg(Yo,Ym,Yc,Ymsk,vx_vol,0,job,seg8t,tis); 
        spm_progress_bar('Set',i - 0.7);
        
  
        % == refine SPM segmentation or just prepare some maps == 
        %  - store the changes also in the output structure
        if job.opts.refine && ~out(i).CTseg
          stime = cat_io_cmd('  Refine SPM','g5','',job.opts.verb>1,stime); 
          [Yc,Ye,Ya,out(i).boney_refine_SPM] = boney_segment_refineSPM(Yo,Ym,Yc,Ya,Ybraindist0,tis,tismri);
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
          boney_segment_extractbone(Vo,Ym,Yc,Ye,Ya,Ymsk,seg8t,tis,out(i),job,vx_vol,YaROIname,RES,BB);
        spm_progress_bar('Set',i - 0.4);
 

        % == restore resolution & boundary box == 
        Ym             = cat_vol_resize( Ym            ,'dereduceV'    ,RES);
        [Ya, Ymsk]     = cat_vol_resize({Ya, Ymsk}     ,'dereduceV'    ,RES,'nearest');
        [Ym, Ymsk, Ya] = cat_vol_resize({Ym, Ymsk, Ya} ,'dereduceBrain',BB); 
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
        %  - although the surfaces can be saved it does not suport surface
        %    registration now!
        if job.opts.bmethod>1
          stime = cat_io_cmd('  Extract bone surfaces','g5','',job.opts.verb>1,stime); 
          [Si, Stm, sROI] = boney_segment_create_bone_surface(Vo, Ybonepp, Ybonemarrow, Ybonethick, Yheadthick, Ya, Ymsk, YaROIname, out(i), job);
        else
          Si = ''; Stm = '';
        end

      else
        % == fast SPM mat-based pipeline ==
        %  - fast version that looks only for (and exports) the bone tissue 
        %    values and creates some report (Yo==Ym?)
        %  - no other tissues (e.g. GM etc.) are loaded and processed!
        [Vo, Ym, Yc, Ybonemarrow, tismri, Si, Stm, Affine] = ...
          boney_segment_fst(P{i}, out(i).P.cls{4}, job, seg8t, tis, vx_vol);
      end
      spm_progress_bar('Set',i - 0.2);
  
     

      % == create output structure
% ############## can some of these fields be removed before?      
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
      out(i).date    = datestr(stime3,'YYYYmmDD-HHMMSS');
      if exist('vROI','var'), out(i).vROI = vROI;   end
      if exist('sROI','var'), out(i).sROI = sROI;   end
      
      
      % == additional combined measures for the bone 
      if exist('vROI','var')
        out(i).main.vBMDH{gi}(si,1)  = -out(si).vROI.bonecortex(3) + out(si).vROI.bonethickness(1) ...
                                       -out(si).classic.bone_med   + out(si).vROI.bonethickness(3); % best for BMDH
      end 
      if exist('sROI','var')
        out(i).main.sBMDH{gi}(si,1)  = -out(si).sROI.bonecortex(3) + out(si).sROI.bonethickness(1) ...
                                       -out(si).classic.bone_med   + out(si).sROI.bonethickness(3);
        out(i).main.sBMDH2{gi}(si,1) = -out(si).sROI.bonecortex(3) + out(si).sROI.bonethickness(3) ...
                                       -out(si).classic.bone_med   + out(si).sROI.bonemarrow(3);
      end       
          


      % export to XML 
      cat_io_xml(out(i).P.xml, out(i)); % export to XML 
      

      %% == create report ==
      if job.output.report
        stime = cat_io_cmd('  Create report','g5','',job.opts.verb>1,stime); 
        boney_segment_print_figure(Vo, Ym, Yc, Ybonemarrow, Si, Stm, out(i), job, Affine);
      end
      if job.opts.verb>1, fprintf('% 5.0fs\n',etime(clock,stime)); end
      spm_progress_bar('Set',i - 0.1);
      rerunstr = '';
      
      
    else
    % == load previous result, reprocess in case of error == 
      try
        outi     = cat_io_xml(out(i).P.xml);
        rerunstr = 'loaded';
      catch 
        optb     = job; optb.verb = 0; optb.rerun = 1; optb.files = job.files(i);
        outi     = boney_segment(optb);
        rerunstr = 'rerun';
      end
      out    = cat_io_mergeStruct(out,outi,[],i); 
    end


    % == prepare output and cleanup files ==
    Pout = boney_segment_cleanup(Pout,out,job,i);


    % == print command line report == 
    boney_segment_cmdline(job,out,i,stime2,rerunstr);
    spm_progress_bar('Set',i);

  end


% ############ final csv-export depending on processing level?
% - use/add flag?
  if job.output.writeCSV && numel(Pout.xml)>1
    %%
    stime = cat_io_cmd('Create CSV with all subjects','g9','',job.opts.verb>0); 
    rlevel = {'boney_default','boney_expert'}; 
    matlabbatch{1}.spm.tools.boney.xml2csv.files        = Pout.xml'; 
    matlabbatch{1}.spm.tools.boney.xml2csv.outdir       = {out(1).P.orgpp};
    matlabbatch{1}.spm.tools.boney.xml2csv.fname        = sprintf('Boney_xml_%s.csv',datestr(stime,'YYYYmmDD-HHMMSS'));
    matlabbatch{1}.spm.tools.boney.xml2csv.fieldnames   = {' '};
    matlabbatch{1}.spm.tools.boney.xml2csv.avoidfields  = {''};
    matlabbatch{1}.spm.tools.boney.xml2csv.report       = rlevel{ ( boned.expertgui > 0) + 1 };

    % run silently 
    txt = evalc('spm_jobman(''run'',matlabbatch);');
    if job.opts.verb>0, fprintf('% 5.0fs\n',etime(clock,stime)); end
  end


  % final print and cleanup
  if job.opts.verb
    spm_progress_bar('Clear');
    fprintf('done.\n')
  end
  spm_progress_bar('Clear'); 
end
