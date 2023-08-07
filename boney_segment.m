function out = boney_segment(job)
%boney_segment. Extract bone measures based on SPM segmentation.
% ... DETAILED DESCRIPTION ...
% 
%  out = boney_segment(job); 
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
%                 - reducion factor vs. vertices in humans:   
%                     1~120k, 2~30k, 3~13k, 4~7k, 6~3k, 8~2k
%                 - robust results for 1-4
%    .mask      .. mask problematic regions (default=1)
%   .output     .. structure of output values
%    .report    .. create JPG report file (default=1) 
%    .writevol  .. write volume output (default=0) 
%    .writesurf .. write surface data (default=0)
%    
%   .out
%    .raw       .. original results (e.g. min,max,median,std, ...) for the 
%                  different regions etc. 
%    .spm8     .. SPM data from the Unified Segmentation  
%    .tis      .. global estimated 
%    .vol      .. volume-based global/regional evaluation of the bone/head
%    .surf     .. surface-based global/regional evaluation of the bone/head
%   .main      .. most relevant values (see ...)
% 
% _________________________________________________________________________
%
% Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________
% $Id$

% TODO: 
%  * bone / bone marrow TPM concept
%  * children case concept
%  * skull-stripping / deface defintions 

  global cat_err_res; %#ok<GVMIS> % for CAT error report
  
  if ~exist('job','var'), help boney_segment; return; end
  
  P = job.files; 
  if isempty(P) || isempty(P{1})
    help boney_segment; 
    return; 
  end


  % == defaults paraemter ==
  def.files             = {};       % SPM m-files
  def.opts.verb         = 2;        % display progress (0 - be silent, 1 - one line per subject, 2 - details)
  def.opts.pmethod      = 1;        % preprocessing method: 1-SPM12, 2-CAT12 
  def.opts.bmethod      = 2;        % method: 0 - SPM seg8t eval only, 1 - volume-based, 2 - surface based
  def.opts.prerun       = 0;        % avoid preprocessing in matching files of SPM/CAT are available
  def.opts.bias         = 1;        % strong bias correction
  def.opts.rerun        = 0;        % rerun processing (0 - no load previous result if available, 1 - yes) 
  def.opts.affreg       = 0;        % do affine registration based on new skull-stripping  
  def.opts.reduce       = 4;        % voxel binning for surface creation (higher values create surfaces with less vertices, i.e., details)
                                    %  - reducion factor vs. vertices in humans:   1~120k, 2~30k, 3~13k, 4~7k, 6~3k, 8~2k
                                    %  - robust results for 1-4
  def.opts.refine       = 1;        % refine segmenation 
  def.opts.subdirs      = 1;        % use subdirectories 
  def.opts.snspace      = [80,7,2]; % cmd print format
  def.opts.normCT       = 0; 
  def.opts.Patlas       = fullfile(spm('dir'),'toolbox','cat12','templates_volumes_not_distributed','cat_KADA_bone-regions.nii');
  def.opts.Pmask        = fullfile(spm('dir'),'toolbox','cat12','templates_volumes_not_distributed','cat_KADA_bone-marrowmask.nii');
  def.opts.MarkColor    = cat_io_colormaps('marks+',40); 
  def.opts.reslim       = 1.5;      % general resolution limit 
  def.opts.expert       = 2;        % user level (0 - default, 1 - expert, 2 - developer)
  def.opts.fst          = 1;        % estimate also first version 
  def.output.report     = 1;        % write report
  def.output.writevol   = 1;        % write volume output 
  def.output.writesurf  = 1;        % write surface data 
  job                   = cat_io_checkinopt(job,def);
  job.output.report     = min(job.output.report,max(1,job.opts.bmethod)); % no surface output without surface processing


  % filenames for dependencies 
  [out,job.opts.fmethod] = boney_segment_filenames(P,job);
  

% #######################
% if ..., return; end % matlabbatch
%
% call SPM / CAT preprocessing
% * use some different parameterization 
% * 
% #######################
  if job.opts.fmethod % if as input-files segments were used we use them and do not reprocess at all 
    boney_segment_preprocessing(P, out, job.opts.pmethod, job.opts.bias, job.opts.prerun);
  end



   
  % =======================================================================
  % TODO: 
  % * Check volumes and CSV to detect bad images (QC?) and segmentations 
  %   (soft values and high variance in values)
  %   >> warning + report + suggestion (affine, resolution)
  %
  % * add real fat measures 
  %   - bone vs. head thickness? no you will need fat vs. skin/musles
  %   - required for better scaling at some point ...
  %
  % * use intnormed image?
  %   - our intensity values are going over different classes - intra vs. intersite variability 
  %   - use CSF-WM contrast (but what about HH?) 
  %
  % * use fat segmentation to avoid critical regions
  %
  % * test SPM progress bar
  % * add rerun parameter
  %
  % * test Brutforce segmentation 1 and 2  
  % * update Yc map (at least print this one) ... diff muscles/bone?
  % * Ym-scaling in MT incrorrect
  % * colorbar vol from reportfig
  % * colorbar histogram for surf
  % * simple version just to report SPM segmentation 
  %   > render skull + brain with pbt 1 mm with thicknes ;-)
  %
  % * list of used CAT functions
  %
  % =======================================================================


  % create table elements for the in-line report
  % TODO: 
  % * Visual feature:  one line progress >> first name than progress 
  % * write table as csv and add conclusion on command line
  [~, Tline, ~, ~ ,MAfn,matm,mmatm] = boney_segment_prepare_print(P,job);
  stime3 = clock; i = 1; %#ok<NASGU> 
  


  %% main loop over all subjects
  spm_progress_bar('Init',numel(job.files),'Bone Segment','Volumes Complete');
  for i = 1:numel(P)
    clearvars -global cat_err_res;
    cat_err_res.stime        = clock; 
    cat_err_res.cat_warnings = cat_io_addwarning('reset'); % reset warnings 
    stime2 = clock;   

    if ... %cat_io_rerun(which(mfilename),out(i).P.xml) || ...
        cat_io_rerun(out(i).P.org,out(i).P.xml,0) || job.output.rerun
  
  
      % == GET SPM DATA ==
      % - get and evaluate SPM preprocessing structure (seg8t vs. tis)
      % - get cmd line output table row (matm & mmatm)
      stime = cat_io_cmd('  Load SPM','g5','',job.opts.verb>1); 
      [ seg8t , tis , matm(i,:), mmatm(i,:), vx_vol ] = boney_segment_get_segmat(out(i),MAfn); 
      job.reslim = max(job.opts.reslim,mean(vx_vol)); % limit processing resolution 

  
      if job.opts.bmethod
        % == load MRI and segmentation ==
        %ds('d2sm','',vx_vol,Ym .* (0.5 + 0.5*Ybone) ,Ya,90)
        stime = cat_io_cmd('  Load MRIs','g5','',job.opts.verb>1,stime); 
        if out(i).CTseg, job.affreg = -1; end 
        [Vo,Yo,Yc,Ya,Ymsk,Ym, Affine, RES, BB] = boney_segment_loadMRI(out(i).P, job, seg8t, tis, 1:5, 25); % CT rand issue
        spm_progress_bar('Set',i - 0.8); vx_vol = RES.vx_volr; 
      
        if job.opts.fst
          stime = cat_io_cmd('  Fast bone measures','g5','',job.opts.verb>1,stime); 
          out(i).fst = boney_segment_simpleBone(seg8t,Yo,Yc); 
        end
  
        % == evaluate/test the SPM segmentation == 
        stime  = cat_io_cmd('  Evaluate MRIs','g5','',job.opts.verb>1,stime); 
        [tismri,Ybraindist0] = boney_segment_evalSPMseg(Yo,Ym,Yc,Ymsk,vx_vol,0,job,seg8t,tis); % full estimation with refined tissue peaks
        spm_progress_bar('Set',i - 0.7);
        
  
        %% == refine SPM segmentation or just prepare some maps == 
        if job.opts.refine && ~out(i).CTseg
          stime = cat_io_cmd('  Refine SPM','g5','',job.opts.verb>1,stime); 
          [Yc,Ye,clsmod] = boney_segment_refineSPM(Yo,Ym,Yc,Ybraindist0,tis,tismri);
        else
          Ye = cell(0); 
        end
        spm_progress_bar('Set',i - 0.6);
      
  
        % == get bone measures == 
        % Y*   .. bone/head maps for surface mapping 
        % vROI .. extracted global/regional bone/head values 
        stime = cat_io_cmd('  Extract bone measures','g5','',job.opts.verb>1,stime);
        [Ybonepp,Ybonethick,Ybonemarrow,Yhdthick,vROI] = boney_segment_extractbone(Vo,Yo,Ym,Yc,Ye,Ya,Ymsk,seg8t,tis,out(i),job,vx_vol,RES,BB);
        spm_progress_bar('Set',i - 0.4);
 

        % restore resolution & boundary box
        [Yo, Ym]   = cat_vol_resize({Yo, Ym}   ,'dereduceV'  ,RES);
        [Ya, Ymsk] = cat_vol_resize({Ya, Ymsk} ,'dereduceV'  ,RES,'nearest');
        [Yo, Ym, Ya, Ymsk] = cat_vol_resize({Yo, Ym, Ya, Ymsk} ,'dereduceBrain',BB); 
        for ci = 1:numel(Yc)
          Yc{ci} = cat_vol_resize( Yc{ci}  ,'dereduceV'    ,RES);
          Yc{ci} = cat_vol_resize( Yc{ci}  ,'dereduceBrain',BB); 
        end

  
if job.opts.bmethod>1
        %% == surface-bases processing ==
        % S*   .. bone surface with bone/head measures
        % sROI .. extracted global/regional bone/head surface values
        stime = cat_io_cmd('  Extract bone surfaces','g5','',job.opts.verb>1,stime); 
        [Si,St,Stm,Sth,sROI] = boney_segment_create_bone_surface(Vo,Ym,Yc,Ybonepp,Ybonethick,Ybonemarrow,Yhdthick,Ya,Ymsk,out(i),job);


        mn = cat_stat_kmeans(Si.facevertexcdata(St.facevertexcdata>5 & Si.facevertexcdata>0 ),3); % maximum value 
        tismri.iBone = mn(2); %cat_stat_nanmean(Si.facevertexcdata);
        tismri.iBonemn = mn;
        [tismri.iBonemn3,tismri.iBonesd3,tismri.iBonevx3] = ...
          cat_stat_kmeans(Si.facevertexcdata(St.facevertexcdata>5 & Si.facevertexcdata>0 ),3); % maximum value 
        [tismri.iBonemn2,tismri.iBonesd2,tismri.iBonevx2] = ...
          cat_stat_kmeans(Si.facevertexcdata(St.facevertexcdata>5 & Si.facevertexcdata>0 ),2); % maximum value 

% same as ROI
tismri.surfmd  = cat_stat_nanmedian( Si.facevertexcdata(St.facevertexcdata>5 & Si.facevertexcdata>0 ) );
tismri.surfmn  = cat_stat_nanmean(   Si.facevertexcdata(St.facevertexcdata>5 & Si.facevertexcdata>0 ) );
tismri.surfsd  = cat_stat_nanstd(    Si.facevertexcdata(St.facevertexcdata>5 & Si.facevertexcdata>0 ) );
tismri.surfiqr = iqr(                Si.facevertexcdata(St.facevertexcdata>5 & Si.facevertexcdata>0 ) );
        % thickness
        mn = cat_stat_kmeans(        St.facevertexcdata(St.facevertexcdata>5 & Si.facevertexcdata>0 ),3); % median thickness to avoid low and high outliers
        tismri.tBone   = mn(2); %cat_stat_nanmean(St.facevertexcdata);
        tismri.tBonemn = mn; 

        [tismri.volmn3, tismri.volsd3, tismri.volvx3]  = cat_stat_kmeans( Ybonemarrow(Yc{4}(:)>.5) ,3); % median thickness to avoid low and high outliers
        [tismri.volmn2, tismri.volsd2, tismri.volvx2]  = cat_stat_kmeans( Ybonemarrow(Yc{4}(:)>.5) ,2); % median thickness to avoid low and high outliers
tismri.volmd   = cat_stat_nanmedian( Ybonemarrow(Yc{4}(:)>.5) );
tismri.volmn   = cat_stat_nanmean(   Ybonemarrow(Yc{4}(:)>.5) );
tismri.volsd   = cat_stat_nanstd(    Ybonemarrow(Yc{4}(:)>.5) );
tismri.voliqr  = iqr(                Ybonemarrow(Yc{4}(:)>.5) );
% head thickness
%Stm.facevertexcdata = Stm.facevertexcdata ./ sum(tismri.vol(1:3))^(1/3);
tismri.headthickkm2 = cat_stat_kmeans(    Stm.facevertexcdata(Stm.facevertexcdata>1 & Stm.facevertexcdata<30 ),2); 
tismri.headthickkm3 = cat_stat_kmeans(    Stm.facevertexcdata(Stm.facevertexcdata>1 & Stm.facevertexcdata<30 ),3); 
tismri.headthickmd  = cat_stat_nanmedian( Stm.facevertexcdata(St.facevertexcdata>5 & Stm.facevertexcdata>0 ) );
tismri.headthickmn  = cat_stat_nanmean(   Stm.facevertexcdata(St.facevertexcdata>5 & Stm.facevertexcdata>0 ) );
tismri.headthicksd  = cat_stat_nanstd(    Stm.facevertexcdata(St.facevertexcdata>5 & Stm.facevertexcdata>0 ) );
tismri.headthickiqr = iqr(                Stm.facevertexcdata(St.facevertexcdata>5 & Stm.facevertexcdata>0 ) );
else
  St = ''; Si = ''; Stm = '';
end

      else
        %[pp,ff] = spm_fileparts(P{i}); ff = ff(2:end);
        Affine = seg8t.Affine; 
        %Vo = seg8t.image;
        Ym = []; Yc = {}; %Ybonemarrow = []; 
        tismri = struct();
      
        %% get bias corrected original image
        if 1
          Vo = spm_vol(P{i});
          Yo = single(spm_read_vols(Vo));
      
          Pc4  = fullfile(out(i).P.orgpp,sprintf('c%d%s%s',4,out(i).P.orgff,out(i).P.ee));
          Vc4  = spm_vol(Pc4); 
          Yc4  = single(spm_read_vols(Vc4));
          [Ytiv,redR] = cat_vol_resize(Yc4,'reduceV',vx_vol,4,8,'max');
          Ytiv = cat_vol_morph(Ytiv>.5,'ldc',4) & Ytiv<.5; 
          Ytiv = cat_vol_morph(Ytiv>.5,'o',1); 
          Ytiv = cat_vol_resize(Ytiv,'dereduceV',redR);
  
          Ybonemarrow = single(Yo/tis.seg8o(3)) .* (Yc4>.5) / tis.seg8o(3);
          
          % get measures
          tismri.TIV     = cat_stat_nansum(Ytiv(:)>0.5) .* prod(vx_vol) / 1000;
          tismri.den     = [nan nan nan cat_stat_nansum(Yc4(:))     .* prod(vx_vol) / 1000 nan nan]; 
          tismri.vol     = [nan nan nan cat_stat_nansum(Yc4(:)>0.5) .* prod(vx_vol) / 1000 nan nan]; 
          tismri.volr    = tismri.vol ./ tismri.TIV;
          tismri.Tth     = [nan nan nan cat_stat_nansum(Yo(Yc4(:)>0.5)) / tis.seg8o(2) nan nan]; 
          tismri.Tsd     = [nan nan nan cat_stat_nanstd(Yo(Yc4(:)>0.5)) / tis.seg8o(2) nan nan]; 
          mn             = cat_stat_kmeans(Ybonemarrow(Yc4(:)>.5),3); % maximum value 
          tismri.iBone   = mn(2); %cat_stat_nanmean(Si.facevertexcdata);
          tismri.iBonemn = mn;
          [tismri.volmn3,tismri.volsd3,tismri.volvx3]  = cat_stat_kmeans( Ybonemarrow(Yc4(:)>.5) ,3); % median thickness to avoid low and high outliers
          [tismri.volmn2,tismri.volsd2,tismri.volvx2]  = cat_stat_kmeans( Ybonemarrow(Yc4(:)>.5) ,2); % median thickness to avoid low and high outliers
          tismri.volmd  = cat_stat_nanmedian( Ybonemarrow(Yc4(:)>.5));
          tismri.volmn  = cat_stat_nanmean( Ybonemarrow(Yc4(:)>.5) );
          tismri.volsd  = cat_stat_nanstd( Ybonemarrow(Yc4(:)>.5) );
          tismri.voliqr = iqr( Ybonemarrow(Yc4(:)>.5) );
%          nout(si).tis.seg8conr; 
         
          %% write output maps
          if job.output.writevol
            %%
            tdim = seg8t.tpm(1).dim; 
            M0   = seg8t.image.mat;          
            M1   = seg8t.tpm(1).mat;
        
            % affine and rigid parameters for registration 
            % if the rigid output is incorrect but affine is good than the Yy caused the problem (and probably another call of this function) 
            R               = spm_imatrix(seg8t.Affine); R(7:9)=1; R(10:12)=0; R=spm_matrix(R);  
            Mrigid          = M0\inv(R)*M1;                                                          % transformation from subject to registration space (rigid)
            Maffine         = M0\inv(seg8t.Affine)*M1;                                                 % individual to registration space (affine)
            mat0a           = seg8t.Affine\M1;                                                         % mat0 for affine output
            mat0r           = R\M1;                                                                  % mat0 for rigid ouput
   
            % settings for new boundary box of the output images 
            trans.native.Vo = seg8t.image(1); 
            trans.native.Vi = seg8t.image(1);
            trans.affine    = struct('odim',tdim,'mat',M1,'mat0',mat0a,'M',Maffine,'A',seg8t.Affine);  % structure for cat_io_writenii
            trans.rigid     = struct('odim',tdim,'mat',M1,'mat0',mat0r,'M',Mrigid ,'R',R);           % structure for cat_io_writenii
          
            job.output.bonemarrow  = struct('native',1,'warped',0,'dartel',3);
            % midline map also for masking masking
            cat_io_writenii(Vo,Ybonemarrow,'',sprintf('bonemarrow%d_',job.opts.bmethod), ...
              'bone percentage position map','uint16',[0,0.001], ... 
              min([1 0 2],[job.output.bonemarrow.native job.output.bonemarrow.warped job.output.bonemarrow.dartel]),trans);
          end
        end
    
        % no surfaces
        St = ''; Si = ''; Stm = '';
      end
      spm_progress_bar('Set',i - 0.2);
      if out(i).CTseg
        mn = cat_stat_kmeans(Yo(Yc{4}>.3),3); % maximum value 
        matm{i,end-3} = median(Yo(Yc{4}>.3)); % classic
      elseif exist('Yc4','var')
        mn = cat_stat_kmeans(Yo(Yc4>.3) / tis.seg8o(2),3); % maximum value 
        matm{i,end-3} = median(Yo(Yc4>.3) / tis.seg8o(2)); % classic
      else
        mn = cat_stat_kmeans(Yo(Yc{4}>.3) / tis.seg8o(2),3); % maximum value 
        matm{i,end-3} = median(Yo(Yc{4}>.3) / tis.seg8o(2)); % classic
      end
      matm{i,end-2}   = mn(3); %cat_stat_nanmean(Si.facevertexcdata);
      try
        matm{i,end-1}   = tismri.tBone; %cat_stat_nanmean(Si.facevertexcdata);
        matm{i,end}     = tismri.headthickmd; %cat_stat_nanmean(Si.facevertexcdata);
      end
  
      %% == create report ==
      if job.output.report
        stime = cat_io_cmd('  Create report','g5','',job.opts.verb>1,stime); 
        boney_segment_print_figure(Vo,Ym,Yc,Ybonemarrow,Si,St,Stm,seg8t,tis,tismri,out(i).P,job,Affine);
      end
      if job.opts.verb>1, fprintf('% 5.0fs\n',etime(clock,stime)); end
      spm_progress_bar('Set',i - 0.1);
    

      
      %% == create output structure
      if isfield(seg8t.image,'private'), seg8t.image = rmfield(seg8t.image,'private'); end
      if isfield(seg8t.tpm,  'private'), seg8t.tpm   = rmfield(seg8t.tpm  ,'private'); end
      if isfield(seg8t, 'sett'), seg8t  = rmfield(seg8t ,'sett'); end
      out(i).spm8    = seg8t; 
      out(i).tis     = tis; 
      out(i).repstr  = struct('mmatm',mmatm(i,:), 'matm', {matm(i,:)}, 'Theader', {MAfn}); 
      if ~isempty('tismri'), out(i).tismri  = tismri; end
      if exist('sROI','var'), out(i).sROI     = sROI;    end
      if exist('vROI','var'), out(i).vROI     = vROI;    end

      % export to XML 
      outxml = out(i);
      if isfield(outxml.spm8,'tpmA'),   outxml.spm8 = rmfield(outxml.spm8,'tpmA'); end
      cat_io_xml(out(i).P.xml,outxml); 
      rerunstr = '';
    else
    % == load prevoius result, reprocess in case of error == 
      try
        outi     = cat_io_xml(out(i).P.xml);
        rerunstr = 'loaded';
      catch 
        optb     = job; optb.verb = 0; optb.rerun = 1; optb.files = job.files(i);
        outi     = boney_segment(optb);
        
        rerunstr = 'rerun';
      end
      out    = cat_io_mergeStruct(out,outi,[],i); 
%      out(i) = out(end); out(end) = []; 
      try
        matm(i,:)  = out(i).repstr.matm; 
        mmatm(i,:) = out(i).repstr.mmatm;
      catch
        optb       = job; optb.verb = 0; optb.rerun = 1; optb.files = job.files(i);
        outi       = boney_segment(optb);
        out        = cat_io_mergeStruct(out,outi); 
        out(i)     = out(end); out(end) = []; 
        matm(i,:)  = out(i).repstr.matm; 
        mmatm(i,:) = out(i).repstr.mmatm;
      end
    end

    if job.opts.verb
      %% Todo
      %  - Link volume with SEG overlay 
      % [-] Volume with surf overlay + "[MBI3D,TH3D]" Links
      % * create a table header / structure!
      % * save processing output mat/xml
      
      % - link processed file (vol + surf?) wherefore?
      % - colorcoded values by CAT intensity map?
      % * add warning Letter!
      % * tiss error with class 
      % * thin bone thick head error?
      %
      if out(i).CTseg
        %%
        matm{i,end-2} = matm{i,end-2} / 1000;
        matm{i,end-1} = matm{i,end-1} / 1000;
        %%
        cat_io_cprintf( job.MarkColor( min(size(job.MarkColor,1),max(1,floor( (matm{i,end-2} - 1) * ...
          size(job.MarkColor,1)))),:),sprintf( Tline,i,...
            sprintf( sprintf('%%%ds',job.snspace(1)-8) , spm_str_manip(P{i},['a' num2str(job.snspace(1) - 14)])), ...
            ... spm_file( sprintf( sprintf('%%%ds',job.snspace(1)-8) , spm_str_manip(P{i},['a' num2str(job.snspace(1) - 14)])),'link',sprintf('spm_display(''%s'');',P{i})), ...
            matm{i,:},etime(clock,stime2)));
      else    
        cat_io_cprintf( job.opts.MarkColor( min(size(job.opts.MarkColor,1),max(1,floor( matm{i,end-2} * 3 / 9.5 * ...
          size(job.opts.MarkColor,1)))),:),sprintf( Tline,i,...
            sprintf( sprintf('%%%ds',job.opts.snspace(1)-8) , spm_str_manip(P{i},['a' num2str(job.opts.snspace(1) - 14)])), ...
            ... spm_file( sprintf( sprintf('%%%ds',job.snspace(1)-8) , spm_str_manip(P{i},['a' num2str(job.snspace(1) - 14)])),'link',sprintf('spm_display(''%s'');',P{i})), ...
            matm{i,:},etime(clock,stime2)));
      end

      fprintf(' %s',rerunstr);

      warn = cat_io_addwarning(1); if numel(warn)>0, cat_io_cprintf('warn',sprintf(' %dW',numel(warn))); end
      errn = cat_io_addwarning(2); if numel(errn)>0, cat_io_cprintf('err' ,sprintf(' %dE',numel(errn))); end
      fprintf('\n');
    end

    spm_progress_bar('Set',i);

  end
  % final print
  if job.opts.verb
    spm_progress_bar('Clear');
    fprintf('done.\n')
  end
end
