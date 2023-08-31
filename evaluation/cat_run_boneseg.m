function out = cat_run_boneseg(opt)
%cat_run_boneseg. Extract bone measures based on SPM segmentation.
% ... DETAILED DESCRIPTION ...
% 
%  out = cat_run_boneseg(opt); 
%
%  opt      .. structure of input values
%   .files     .. SPM m-files
%   .method    .. method used for evaluation (default=2): 
%                  0 - SPM seg8t mat values evaluation (fast)
%                  1 - SPM volumes with problems for high intensity bone 
%                      marrow 
%                  2 - refined volumes  
%   .writevol  .. write volume output (default=1) 
%   .writesurf .. write surface data (default=1)
%   .report    .. create JPG report file (default=1) 
%   .verb      .. display progress 
%                 (0 - be silent, 1 - one line per subject, 2 - details)
%   .affreg    .. do affine registration based (default=1) 
%   .reduce    .. voxel binning for surface creation (higher values create 
%                 surfaces with less vertices, i.e., details)
%                 - reducion factor vs. vertices in humans:   
%                     1~120k, 2~30k, 3~13k, 4~7k, 6~3k, 8~2k
%                 - robust results for 1-4
%   .mask      .. mask problematic regions (default=1)
%
%  out         .. structure of output values
%   .raw       .. original results (e.g. min,max,median,std, ...) for the 
%                 different regions etc. 
%    .spm8     .. SPM data from the Unified Segmentation  
%    .tis      .. global estimated 
%    .vol      .. volume-based global/regional evaluation of the bone/head
%    .surf     .. surface-based global/regional evaluation of the bone/head
%   .main      .. most relevant values (see ...)
% 
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
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
  
  if ~exist('opt','var'), return; end
  
  P = opt.files; 
  if isempty(P) || isempty(P{1})
    help cat_run_bonseg; 
    return; 
  end


  % == defaults ==
  def.files     = {};       % SPM m-files
  def.verb      = 2;        % display progress (0 - be silent, 1 - one line per subject, 2 - details)
  def.method    = 2;        % method: 0 - SPM seg8t eval only, 1 - SPM vols, 2 - refined SPM vols,  
  def.report    = 1;        % write report
  def.rerun     = 0;        % rerun processing (0 - no load previous result if available, 1 - yes) 
  def.subdirs   = 0;        % use subdirs 
  def.writevol  = 1;        % write volume output 
  def.writesurf = 1;        % write surface data 
  def.affreg    = 0;        % do affine registration based on new skull-stripping  
  def.reduce    = 4;        % voxel binning for surface creation (higher values create surfaces with less vertices, i.e., details)
                            %  - reducion factor vs. vertices in humans:   1~120k, 2~30k, 3~13k, 4~7k, 6~3k, 8~2k
                            %  - robust results for 1-4
  def.reslim    = 1.5;      % general resolution limit 
  def.mask      = 1;        % mask problematic regions
  def.snspace   = [80,7,2]; % cmd print format
  def.expert    = 2;        % user level (0 - default, 1 - expert, 2 - developer)
  def.fst       = 1;        % estimate also first version 
  def.normCT    = 0; 
  def.Patlas    = fullfile(spm('dir'),'toolbox','cat12','templates_volumes_not_distributed','cat_KADA_bone-regions.nii');
  def.Pmask     = fullfile(spm('dir'),'toolbox','cat12','templates_volumes_not_distributed','cat_KADA_bone-marrowmask.nii');
  %def.Pmask     = fullfile(spm('dir'),'tpm','labels_Neuromorphometrics_bones.nii');
  def.MarkColor = cat_io_colormaps('marks+',40); 
  opt           = cat_io_checkinopt(opt,def);
  
  % filenames for dependencies 
  out = bonefilenames(P,opt);
  % if ..., return; end % matlabbatch

    
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
  
  
  
  % create table elements for the in-line report
  % TODO: 
  % * one line >> first name than progress 
  % * better integration of links and seperation characters?
  % * Pll - not really interesting 
  % * cls numbers? - most time snot relevant > but for report
  % * write table as csv and add conclusion on command line
  [Theader, Tline,Tavg, ~ ,MAfn,matm,mmatm] = prepare_print(P,opt);
  stime3 = clock; i = 1; 
  
  %%
  spm_progress_bar('Init',numel(opt.files),'Bone Segment','Volumes Complete');
  for i = 1:numel(P)
    clearvars -global cat_err_res;
    cat_err_res.stime        = clock; 
    cat_err_res.cat_warnings = cat_io_addwarning('reset'); % reset warnings 
    stime2 = clock;   

    if ... %cat_io_rerun(which(mfilename),out(i).P.xml) || ...
        cat_io_rerun(out(i).P.org,out(i).P.xml,0) || opt.rerun
  
  
      % == GET SPM DATA ==
      % - get and evaluate SPM preprocessing structure (seg8t vs. tis)
      % - get cmd line output table row (matm & mmatm)
      stime = cat_io_cmd('  Load SPM','g5','',opt.verb>1); 
      [ seg8t , tis , matm(i,:), mmatm(i,:), vx_vol ] = get_segmat(out(i),MAfn); 
      opt.reslim = max(opt.reslim,mean(vx_vol)); % limit processing resolution 

  
      if opt.method
        % == load MRI and segmentation ==
        %ds('d2sm','',vx_vol,Ym .* (0.5 + 0.5*Ybone) ,Ya,90)
        stime = cat_io_cmd('  Load MRIs','g5','',opt.verb>1,stime); 
        if out(i).CTseg, opt.affreg = -1; end 
        [Vo,Yo,Yc,Ya,Ymsk,Ym, Affine, RES, BB] = loadMRI(out(i).P, opt, seg8t, tis, 1:5, 25); % CT rand issue
        spm_progress_bar('Set',i - 0.8); vx_vol = RES.vx_volr; 
      
        if opt.fst
          stime = cat_io_cmd('  fst bone measures','g5','',opt.verb>1,stime); 
          out(i).fst = simpleBone(seg8t,Yo,Yc); 
        end
  
        % == evaluate/test the SPM segmentation == 
        stime  = cat_io_cmd('  Evaluate MRIs','g5','',opt.verb>1,stime); 
        [tismri,Ybraindist0] = evalSPMseg(Yo,Ym,Yc,Ymsk,vx_vol,0,opt,seg8t,tis); % full estimation with refined tissue peaks
        spm_progress_bar('Set',i - 0.7);
        
  
        %% == refine SPM segmentation or just prepare some maps == 
        if opt.method > 1  && ~out(i).CTseg
          stime = cat_io_cmd('  Refine SPM','g5','',opt.verb>1,stime); 
          [Yc,Ye,clsmod] = refineSPM(Yo,Ym,Yc,Ybraindist0,tis,tismri);
        else
          Ye = cell(0); 
        end
        spm_progress_bar('Set',i - 0.6);
      
  
        % == get bone measures == 
        % Y*   .. bone/head maps for surface mapping 
        % vROI .. extracted global/regional bone/head values 
        stime = cat_io_cmd('  Extract bone measures','g5','',opt.verb>1,stime);
        [Ybonepp,Ybonethick,Ybonemarrow,Yhdthick,vROI] = extractbone(Vo,Yo,Ym,Yc,Ye,Ya,Ymsk,seg8t,tis,out(i),opt,vx_vol,RES,BB);
        spm_progress_bar('Set',i - 0.4);
 

        % restore resolution & boundary box
        [Yo, Ym]   = cat_vol_resize({Yo, Ym}   ,'dereduceV'  ,RES);
        [Ya, Ymsk] = cat_vol_resize({Ya, Ymsk} ,'dereduceV'  ,RES,'nearest');
        [Yo, Ym, Ya, Ymsk] = cat_vol_resize({Yo, Ym, Ya, Ymsk} ,'dereduceBrain',BB); 
        for ci = 1:numel(Yc)
          Yc{ci} = cat_vol_resize( Yc{ci}  ,'dereduceV'    ,RES);
          Yc{ci} = cat_vol_resize( Yc{ci}  ,'dereduceBrain',BB); 
        end

  
        %% == surface-bases processing ==
        % S*   .. bone surface with bone/head measures
        % sROI .. extracted global/regional bone/head surface values
        stime = cat_io_cmd('  Extract bone surfaces','g5','',opt.verb>1,stime); 
        [Si,St,Stm,Sth,sROI] = create_bone_surface(Vo,Ym,Yc,Ybonepp,Ybonethick,Ybonemarrow,Yhdthick,Ya,Ymsk,out(i),opt);



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
        %[pp,ff] = spm_fileparts(P{i}); ff = ff(2:end);
        Affine = seg8t.Affine; 
        %Vo = seg8t.image;
        Ym = []; Yc = {}; %Ybonemarrow = []; 
        tismri = struct();
      
        %% get bias corrected original image
        if 1
          Vo = spm_vol(P{i});
          Yo = single(spm_read_vols(Vo));
      
          [pp,ff,ee] = spm_fileparts(P{i});
          Pc4  = fullfile(pp,sprintf('c%d%s%s',4,ff(2:end),ee));
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
          if opt.writevol
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
            cat_io_writenii(Vo,Ybonemarrow,'',sprintf('bonemarrow%d_',opt.method), ...
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
      matm{i,end-1}   = tismri.tBone; %cat_stat_nanmean(Si.facevertexcdata);
      matm{i,end}     = tismri.headthickmd; %cat_stat_nanmean(Si.facevertexcdata);
    
  
      %% == create report ==
      if opt.report
        stime = cat_io_cmd('  Create report','g5','',opt.verb>1,stime); 
        print_figure(Vo,Ym,Yc,Ybonemarrow,Si,St,Stm,seg8t,tis,tismri,out(i).P,opt,Affine);
      end
      if opt.verb>1, fprintf('% 5.0fs\n',etime(clock,stime)); end
      spm_progress_bar('Set',i - 0.1);
    

      
      %% == create output structure
      if isfield(seg8t.image,'private'), seg8t.image = rmfield(seg8t.image,'private'); end
      if isfield(seg8t.tpm,  'private'), seg8t.tpm   = rmfield(seg8t.tpm  ,'private'); end
      if isfield(seg8t, 'sett'), seg8t  = rmfield(seg8t ,'sett'); end
      out(i).spm8    = seg8t; 
      out(i).tis     = tis; 
      out(i).repstr  = struct('mmatm',mmatm(i,:), 'matm', {matm(i,:)}, 'Theader', {MAfn}); 
      if ~isempty('tismri'), out(i).tismri  = tismri; end
      if exist('vROI','var'), out(i).sROI     = sROI;    end
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
        optb     = opt; optb.verb = 0; optb.rerun = 1; optb.files = opt.files(i);
        outi     = cat_run_boneseg(optb);
        
        rerunstr = 'rerun';
      end
      out    = cat_io_mergeStruct(out,outi,[],i); 
%      out(i) = out(end); out(end) = []; 
      try
        matm(i,:)  = out(i).repstr.matm; 
        mmatm(i,:) = out(i).repstr.mmatm;
      catch
        optb       = opt; optb.verb = 0; optb.rerun = 1; optb.files = opt.files(i);
        outi       = cat_run_boneseg(optb);
        out        = cat_io_mergeStruct(out,outi); 
        out(i)     = out(end); out(end) = []; 
        matm(i,:)  = out(i).repstr.matm; 
        mmatm(i,:) = out(i).repstr.mmatm;
      end
    end

    if opt.verb
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
        cat_io_cprintf( opt.MarkColor( min(size(opt.MarkColor,1),max(1,floor( (matm{i,end-2} - 1) * ...
          size(opt.MarkColor,1)))),:),sprintf( Tline,i,...
            sprintf( sprintf('%%%ds',opt.snspace(1)-8) , spm_str_manip(P{i},['a' num2str(opt.snspace(1) - 14)])), ...
            ... spm_file( sprintf( sprintf('%%%ds',opt.snspace(1)-8) , spm_str_manip(P{i},['a' num2str(opt.snspace(1) - 14)])),'link',sprintf('spm_display(''%s'');',P{i})), ...
            matm{i,:},etime(clock,stime2)));
      else    
        cat_io_cprintf( opt.MarkColor( min(size(opt.MarkColor,1),max(1,floor( matm{i,end-2} * 3 / 9.5 * ...
          size(opt.MarkColor,1)))),:),sprintf( Tline,i,...
            sprintf( sprintf('%%%ds',opt.snspace(1)-8) , spm_str_manip(P{i},['a' num2str(opt.snspace(1) - 14)])), ...
            ... spm_file( sprintf( sprintf('%%%ds',opt.snspace(1)-8) , spm_str_manip(P{i},['a' num2str(opt.snspace(1) - 14)])),'link',sprintf('spm_display(''%s'');',P{i})), ...
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
  if opt.verb
    spm_progress_bar('Clear');
    fprintf('done.\n')
  end
end
%=======================================================================
function out = bonefilenames(P,opt)                                        % ok 
%bonefilenames(P). Pepare output filesnames. 

  for i = 1:numel(P)

    % use subdirectories 
    if opt.subdirs
      out(i).P.mridir    = 'mri'; 
      out(i).P.surfdir   = 'surf';
      out(i).P.reportdir = 'report';
    else
      out(i).P.mridir    = ''; 
      out(i).P.surfdir   = ''; 
      out(i).P.reportdir = ''; 
    end
  
    
    % get original file name based on the type of intput segmentation
    [pp,ff,ee]      = spm_fileparts(P{i}); 
    out(i).CTseg    = contains(ff,'_CTseg');
    if out(i).CTseg
      if strcmp(ff(1:2),'c0'), ffs = 18; elseif strcmp(ff(1:2),'ss'), ffs = 4; end 
      ffe = numel(ff) - 6;
    else
      if ff(1)=='m', ffs = 2; elseif ff(1)=='c', ffs = 3; end 
      ffe = numel(ff);
    end
    out(i).P.orgpp  = pp; 
    out(i).P.orgff  = ff(ffs:ffe);
    out(i).P.ppff   = ff;
    out(i).P.ee     = ee; 

    % input volumes
    out(i).P.org    = fullfile(pp,[ff(ffs:ffe) ee]);
    if out(i).CTseg
      out(i).P.bc   = out(i).P.org;
      out(i).P.seg8 = fullfile(pp,sprintf('mb_fit_CTseg.mat'));
    else
      out(i).P.bc   = P{i};
      out(i).P.seg8 = fullfile(pp,sprintf('%s_seg8.mat',out(i).P.orgff));
    end
    if ~exist(out(i).P.seg8,'file')
      cat_io_cprintf('err','Cannot process "%s" because the seg8.mat is missing. \n',P{i});
      out(i).process = 0;
    end

    % output dirs
    out(i).P.mripath    = fullfile(pp,out(i).P.mridir); 
    out(i).P.surfpath   = fullfile(pp,out(i).P.surfdir); 
    out(i).P.reportpath = fullfile(pp,out(i).P.reportdir); 
    % create dirs if required
    if ~exist(out(i).P.mripath   ,'dir'), mkdir(out(i).P.mripath); end
    if ~exist(out(i).P.surfpath  ,'dir'), mkdir(out(i).P.surfpath); end
    if ~exist(out(i).P.reportpath,'dir'), mkdir(out(i).P.reportpath); end


    % xml/mat output
    out(i).P.report = fullfile(out(i).P.reportpath,sprintf('bonereport%d_%s.jpg',opt.method,ff(2:end)));              
    out(i).P.xml    = fullfile(out(i).P.reportpath,sprintf('catbones%d_%s.xml',opt.method,ff(2:end)));
    out(i).P.mat    = fullfile(out(i).P.reportpath,sprintf('catbones%d_%s.mat',opt.method,ff(2:end)));
    

    % vols
    if opt.writevol
      vols = {'pp','bonemarrow','headthick'};
      for vi = 1:numel(vols)
        % ... subject affine warped ... 
        prefix = {'','r'}; postfix = {'','_affine'};
        for pi = 1:numel(prefix)
          out(i).P.([prefix{pi} vols{vi} postfix{pi}]) = fullfile(out(i).P.mripath, ...
            sprintf('%s%s%d_%s%s%s', prefix{pi}, vols{vi}, opt.method, ff(2:end), postfix{pi}, ee));
        end
      end
    end

    
    % surfs
    if opt.writesurf
      out(i).P.central   = fullfile(out(i).P.surfpath,sprintf('%s%d.central.%s.gii','bone',opt.method,ff));        % central
      out(i).P.marrow    = fullfile(out(i).P.surfpath,sprintf('%s%d.marrow.%s'     ,'bone',opt.method,ff));        % marrow
      out(i).P.thick     = fullfile(out(i).P.surfpath,sprintf('%s%d.thickness.%s'  ,'bone',opt.method,ff));        % thickness
      out(i).P.headthick = fullfile(out(i).P.surfpath,sprintf('%s%d.thickness.%s'  ,'head',opt.method,ff));        % thickness
      out(i).P.marrowmin = fullfile(out(i).P.surfpath,sprintf('%s%d.marrowmin.%s'  ,'bone',opt.method,ff));        % minimal bone values
      out(i).P.marrowmax = fullfile(out(i).P.surfpath,sprintf('%s%d.marrowmax.%s'  ,'bone',opt.method,ff));        % maximum bone values
    end
  end
end
%=======================================================================
function print_figure(Vo,Ym,Yc,Ybonemarrow,Si,St,Stm,seg8t,tis,tismri,Po,opt,Affine) % many things 
%print_figure. Print final report figure with volumes and surfaces.

% ToDo
% * V1: nearest interpolation with optimized data ?
% * V2: low fat issue 
% * V2: pure WM segments but jet for bone?
% * V2: countours?
% 
% * S: altas rendering (only highres) > option?
% * S: atlas nan filtering 
%

%#ok<*AGROW> 

  % in case of updates
  seg8t.Affine = Affine; 
  vx_vol = tis.res_vx_vol;
  clear -global st
  global st %#ok<GVMIS> 

  % fontsettings
  crange    = 8; 
  fontsize  = 10; 
  fontcolor = [0 0 0];
  fontname  = 'monospace';

  % colormaps
  clscol = [ lines(7) ; gray(5)];
  clscol = clscol([5,2,1,4,3,10],:); 
  grycol = gray(10);
  tiscol = cat_io_colormaps('BCGWHw',200);
  jetcol = jet(100)/1.5;
  mrkcol = cat_io_colormaps('marks',100)/1.5;
  
  % setup SPM figure
  fg = spm_figure('FindWin','Graphics'); 
  set(0,'CurrentFigure',fg)
  spm_figure('Clear',fg);
  if isempty(fg)
    fg = spm_figure('Create','Graphics','visible','on'); 
  end
  colormap gray; 
  

  % main text report box with header (filename)
  ax = axes('Position',[0.01 0.75 0.98 0.245],'Visible','off','Parent',fg);
  text(0,0.99, ['Bone marrow extraction: ' ...
    strrep( strrep( spm_str_manip(Vo.fname,'k80d'),'\','\\'), '_','\_') '       '],...
    'FontSize',fontsize+1,'FontWeight','Bold','Interpreter','tex','Parent',ax);




  % == table 1: full spm classes values - maybe only expert? ==
  notab1 = 1; 
  if opt.expert > 1 && ~notab1
    mgid = find(seg8t.mg > 2/numel(seg8t.lkp));
    FN  = {'lkp','mg','mn','vr'};
    FNn = {sprintf('%d/%d SPM-Classes ',numel(mgid),numel(seg8t.lkp)),'Propotion','Mean','Std'}; 
    FNt = {'%d','%0.2f','%0.2f','%0.2f'};
    % gray table rows
    for fni = 1:numel(FN)
      if fni==1
        text(0.01,0.95-            0.035,repmat(' ',1,1205), ...
          'BackgroundColor',[0.94 0.94 .94],'Parent',ax,'FontSize',fontsize*.6); 
      elseif mod(fni,2)
        text(0.01,0.95-(0.035*fni)-0.035,repmat(' ',1,1205), ...
          'BackgroundColor',[0.96 0.96 .96],'Parent',ax,'FontSize',fontsize*.6); 
      end
    end
    % add table content
    for fni = 1:numel(FN)
      seg8text(fni,1)   = text(0.01, 0.96-(0.045*fni) +0.01*(fni==1), ...
        sprintf('\\bf%s',FNn{fni}) ,'Fontname',fontname,'FontSize',...
        fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax); 
      if fni==1, seg8text(fni,1).FontWeight = 'bold'; norm = 1; else, norm = seg8t.(FN{fni})(2); end
      for lkpi = 1:numel(mgid) %numel(seg8t.lkp)
        seg8text(fni,mgid(lkpi)+1) = text(0.11 + 0.87/numel(mgid) * (lkpi), 0.95-(0.045*fni) +0.01*(fni==1), sprintf( FNt{fni}, seg8t.(FN{fni})(mgid(lkpi))/norm ),...
          'Fontname',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax,'HorizontalAlignment','right');
        if fni==1 
          seg8text(fni,mgid(lkpi)+1).FontWeight = 'bold'; seg8text(fni,mgid(lkpi)+1).Color = clscol(seg8t.lkp(mgid(lkpi)),:);  
          seg8text(fni,mgid(lkpi)+1).String = [seg8text(fni,mgid(lkpi)+1).String sprintf('(%d)',sum(seg8t.lkp == seg8t.lkp(mgid(lkpi))) )];
        end
        if fni==2, seg8text(fni,mgid(lkpi)+1).Color = grycol( round(8 - 7*seg8t.mg(mgid(lkpi))) ,:);  end
        if fni==3 
          try
            col = tiscol( max(0, min(200,round(100*seg8t.mn(mgid(lkpi))/max(seg8t.mn(2))))) ,:); 
            seg8text(fni,mgid(lkpi)+1).Color = max(0,min(1,col - max(.3,min(sum(col,2)-2)/4)));  
          end
        end
        if fni==4, seg8text(fni,mgid(lkpi)+1).Color = grycol( min(10,max(1,round(8 - 7 * seg8t.vr(mgid(lkpi)) / max(seg8t.vr(1:3))))) ,:);  end
      end
    end

    table1offset = 0.27;
  else 
    table1offset = 0.02; 
  end



  % == table 2: spm classes values ==
  tis.clsn = {'1-GM','2-WM','3-CSF','4-bone','5-head','6-BG'};
  for ci = 1:min(6,max(seg8t.lkp))
    % indicate cases were a Gaussian is used onyl for a small volume, what
    % is pointing to a to large number of Gaussian
    % ... would have to consider also the difference between Gaussians
    if sum(seg8t.lkp==ci) > max(1,sum(seg8t.mg(seg8t.lkp==ci)'<.2 | abs(gradient(seg8t.mn(seg8t.lkp==ci) ./ tis.WMth)) < .2)) 
      mark = '\color[rgb]{0,0.5,1}';
    else
      mark = ''; 
    end
    tis.clsG{ci} = sprintf('%d %s(%d)', sum(seg8t.lkp==ci), mark, ...
      max(1,sum(seg8t.mg(seg8t.lkp==ci)'<.2 | abs(gradient(seg8t.mn(seg8t.lkp==ci) ./ tis.WMth)) < .2)) );  
  end
  % #####################
  % * add QC value for tissue ? 
  % *
  % #####################
  % density ~ volume and so I use only one measure
  if seg8t.isCTseg && ~opt.normCT 
    tis.seg8x = tismri.Tth; 
    FN  = {'clsn','seg8x','vol','volr'};                                  % 'den','clsG','seg8nv',
    FNn = {'TPM class','Med.Int.','Volume mm3','Volume (%TIV)'};     % 'Volume density','n-Class(well)','Main-std',
    FNt = {'%s','%0.0f','%0.0f','%0.1f'};                                 % '%0.0f','%s','%0.3f',
  else
    FN  = {'clsn','seg8n','vol','volr'};                                  % 'den','clsG','seg8nv',
    FNn = {'TPM class','Norm.Med.Int.','Volume mm3','Volume (%TIV)'};     % 'Volume density','n-Class(well)','Main-std',
    FNt = {'%s','%0.3f','%0.0f','%0.1f'};                                 % '%0.0f','%s','%0.3f',
  end
  if isempty(tismri) || ~isfield(tismri,'vol')
    FN = FN(1:3); FNn = FNn(1:3); FNt = FNt(1:3);
  else
    tis.vol  = tismri.vol; 
    tis.den  = tismri.den;
    tis.volr = tismri.vol ./ sum(tismri.vol(1:3)) * 100; 
    % do not show BG volumes 
    tis.vol(6)  = nan; 
    tis.volr(6) = nan;
  end
  for fni = 1:numel(FN)
    if fni==1
      text(0.01,0.95-table1offset-            0.035,repmat(' ',1,1000 + 205),...
        'BackgroundColor',[0.94 0.94 .94],'Parent',ax,'FontSize',fontsize*.6); 
    elseif mod(fni,2)
      text(0.01,0.95-table1offset-(0.045*fni)-0.0,repmat(' ',1,1000 + round(205*(.6/.4))),...
        'BackgroundColor',[0.96 0.96 .96],'Parent',ax,'FontSize',fontsize*.4); 
    end
  end
  for fni = 1:numel(FN)
    tistext(fni,1)   = text(0.01              , 0.95-table1offset-(0.045*fni) +0.01*(fni==1), ...
      sprintf('\\bf%s',FNn{fni}) ,'Fontname',fontname,'FontSize',...
      fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax);
    if fni==1, tistext(fni,1).FontWeight = 'bold'; end
    for lkpi=numel(tis.clsn):-1:1 % down to avoid overlap
      if iscell(tis.(FN{fni})), tmpval = tis.(FN{fni}){lkpi}; else, tmpval = tis.(FN{fni})(lkpi); end
      if ~isnan(tmpval)
        tistext(fni,lkpi+1) = text(0.12 + 0.082*(lkpi), 0.95-table1offset-(0.045*fni) +0.01*(fni==1), sprintf( FNt{fni}, tmpval ),...
          'Fontname',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax,'HorizontalAlignment','right');
        if fni==1, tistext(fni,lkpi+1).FontWeight = 'bold'; tistext(fni,lkpi+1).Color = clscol(lkpi,:);  end
      end
      %if fni==3, tistext(fni,lkpi+1).Color = grycol( max(1,min(10,round(8 - tmpval))) ,:);  end
    end
    % ###################
    % * add colorrating for tissues 
    % * density vs. volume => error?
    % ###########
  end



  % == table 3: processing and qc parameter? 
  % ... processing options (method-SPM,Bone,Surf)
  % - processing [SPM,CTseg]
  % - number spm-peaks [1 1 2 3 4 2]
  % - method [mat,SPM,ref]
  % 

  % == table 4: imaging and bone parameters
  highBGs = {'low','high','high(MT)','mid'}; %twcol = [0 0 0; 0 0 0.5; 0 .5 0; .5 0 0]; 
  tis.highBGn = highBGs{tis.highBG +1}; % ##########################
  tis.res_RES = mean(vx_vol.^2 )^.5;
  if ~isempty(St)
    tis.iBone   = cat_stat_nanmean(Si.facevertexcdata(Si.facevertexcdata<100));
    tis.thBone  = cat_stat_nanmean(St.facevertexcdata(St.facevertexcdata<100));
  end
  tis.medBone0 = 0; 
  try,  tis.headthickmn = tismri.headthickmn; end
% ################## 
%  * add skull-Stripping and defacing value later 
%  * better definition of low fat intensity rating 
%  * do QC with CNR rating? 
% ################################# 
  FN{1}  = {'weightingn', 'res_RES' , 'highBGn' , 'headFatTypen', 'boneIntTypen'}; 
  FNn{1} = {'Tw'        , 'RES'     , 'BGtype'  , 'fatint',       'fatbone'     };
  FNt{1} = {'%s'        , '%0.2f'   , '%s'      , '%s',           '%s'          }; 

  if seg8t.isCTseg
    tis.minBone = tismri.iBonemn3(1);
    tis.medBone = tismri.iBonemn3(2);
    tis.maxBone = tismri.iBonemn3(3);
    FN{2}  = {'minBone'   , 'medBone' , 'maxBone' , 'thBone'  , 'headthickmn'};
    FNn{2} = {'lBone*'    , 'mBone*'  , 'hBone*'  , 'BthBone' , 'thhead'};
    FNt{2} = {'%0.0f'     , '%0.0f'   , '%0.0f'   , '%0.1f mm', '%0.1f mm'};
  else
    FN{2}  = {'minBone'   , 'medBone' , 'maxBone' , 'medBone0', 'iBone'       , 'thBone'  , 'headthickmn'};
    FNn{2} = {'lBone'     , 'mBone'   , 'hBone'   , 'medBone0', 'iBone'       , 'BthBone' , 'thhead'};
    FNt{2} = {'%0.3f'     , '%0.3f'   , '%0.3f'   , '%0.3f'   , '%0.3f'       , '%0.1f mm', '%0.1f mm'};
  end
  if isempty(St)
    FN{2}(end-1:end) = []; FNn{2}(end-1:end) = []; FNt{2}(end-1:end) = []; 
  end
  for fnj = 1:2
    text(0.01,0.64 - (fnj-1)*0.12 - table1offset,repmat(' ',1,1000 + 120),'BackgroundColor',[0.94 0.94 .94],'Parent',ax,'FontSize',fontsize*.6); 
    for fni = numel(FN{fnj}):-1:1
      segtext(fni,lkpi)   = text(0.03 + 0.075*(fni-1), 0.64-(fnj-1)*0.12-table1offset-(0.00) , ...
        sprintf('\\bf%s',FNn{fnj}{fni}) ,'Fontname',fontname,'FontSize',...
        fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax,'FontWeight','bold','HorizontalAlignment','center');
      if ~isempty( FNt{fnj}{fni} )
        segtext(fni,lkpi+1) = text(0.03 + 0.075*(fni-1), 0.64-(fnj-1)*0.12-table1offset-(0.05), sprintf( FNt{fnj}{fni},  tis.(FN{fnj}{fni}) ),...
          'Fontname',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','tex','HorizontalAlignment','center','Parent',ax);
      end
%      if fni==1, segtext(fni,lkpi+1).Color = twcol(  max(1,min(size(twcol,1),tis.(FN{fnj}{fni}(1:end-1)))) , :); end
      if fni==2, segtext(fni,lkpi+1).Color = mrkcol( max(1,min(100,round( 100 * (tis.(FN{fnj}{fni})+1) / 4 ))) , :); end % RES
      if fnj == 1
        if fni >2 && fni < 6
          switch segtext(fni,lkpi+1).String
            case 'low',  segtext(fni,lkpi+1).Color = [ .0 .0 1.0]; 
            case 'mid',  segtext(fni,lkpi+1).Color = [ .6 .5  .0]; 
            case 'high', segtext(fni,lkpi+1).Color = [ .8 .2  .0]; 
            otherwise,   segtext(fni,lkpi+1).Color = [1.0 .0  .0]; 
          end
        end
      else
        if seg8t.isCTseg
          if fni > 0 && fni<4, segtext(fni,lkpi+1).Color = jetcol( max(1,min(100,round(0.07 * tis.(FN{fnj}{fni})))),:);  end
          if fni > 3,          segtext(fni,lkpi+1).Color = jetcol( max(1,min(100,round(5 * tis.(FN{fnj}{fni})))),:);  end
        else
          if fni > 4,          segtext(fni,lkpi+1).Color = jetcol( max(1,min(100,round(5 * tis.(FN{fnj}{fni})))),:);  end
        end
      end
    end
  end



  %% == error messages 
  % ... show only in case of no refinement or as note?
  if 0 
    wstr   = ''; % empty line
    wname  = {'note','warning','alert'}; 
    wcolor = [0 0 1; 1 0.5 0; 0.8 0 0]; 
    for wi=2:-1:0
      warn = cat_io_addwarning(wi); 
      wn   = numel(warn);
      if wn
        if cat_get_defaults('extopts.expertgui') > -wi
          msg  =  strrep(warn(1).identifier,'_','\_') ; 
          for wmi=2:wn
            msg = [msg ', ' strrep( warn(wmi).identifier,'_','\_') ];
          end
          
          if wn~=1, wnamepl='s'; else, wnamepl=''; end
          wstr = [wstr ...
            sprintf('\\color[rgb]{%0.1f %0.1f %0.1f}\\bf{%d %s%s} (%s)',...
              wcolor(wi+1,:), wn,wname{wi+1},wnamepl, msg)]; 
        end
      end
    end
    text(0.01,0.1,wstr,'Interpreter','tex'); 
  end
 

  % method 0 only read the header and we can print only the basic values/tables 
  if 0 % opt.method == 0
    %%    print(fg, '-djpeg', Po.report); 
    figfs10 = [ findobj(fg,'FontSize',10); findobj(fg,'FontSize',9); findobj(fg,'FontSize',11); findobj(fg,'FontSize',8)]; 
    for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize * .8; end; end %#ok<TRYNC> 
    saveas(fg,Po.report,'png')
    for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize / .8; end; end %#ok<TRYNC> 
    return; 
  end


  % histogram 
  if opt.report > 1 && opt.method > 0 
    %%
    ax2 = axes('Position',[0.61 0.855-table1offset/3 0.40 0.13],'Parent',fg,'Color',[1 1 1]); hold on;
    ax2.YColor = [1 1 1]; ax2.XColor = [1 1 1]; ax2.YTick = []; ax2.XTick = []; 
    ax2 = axes('Position',[0.62 0.855-table1offset/3 0.37 0.13],'Parent',fg,'Color',[1 1 1]); hold on;
    if seg8t.isCTseg && ~opt.normCT 
      for ci = numel(Yc)-1:-1:1
        % (tismri.Tth(2) ./ tis.WMth) correction to print the SPM 
        hhst(ci) = histogram(ax2,Ym( Yc{ci}(:)>.5 & Ym(:)>-1000 & Ym(:)<2000) ,-1000:10:2000, ...
          'LineStyle','none','FaceColor',clscol(ci,:));
        hstmax(ci) = max( hhst(ci).Values )  / 3;
      end
      xlim([-300 1800]); %title('normalized intensities');
    elseif opt.method > 0 
      for ci = numel(Yc)-1:-1:1
        % (tismri.Tth(2) ./ tis.WMth) correction to print the SPM 
        hhst(ci) = histogram(ax2,Ym( Yc{ci}(:)>.5 & Ym(:)>0 & Ym(:)<2) ,0:0.01:2, ...
          'LineStyle','none','FaceColor',clscol(ci,:));
        hstmax(ci) = max( hhst(ci).Values ) ;
      end
      xlim([0 2]); ax2.XTick = 0:1/6:2; % grid on;  %title('normalized intensities');
    else
      % not loaded
      hhst(4) = histogram(ax2,Ybonemarrow(Ybonemarrow(:)>0),0:0.01:2, ...
        'LineStyle','none','FaceColor',clscol(ci,:));
      hstmax(4) = max( hhst(4).Values ) ;
      xlim([0 2]); ax2.XTick = 0:1/6:2; % grid on;  %title('normalized intensities');
    end
    %%
    ylim([0 max(hstmax(1:5)*1.3 )])
    box on; 
    if ~(seg8t.isCTseg && ~opt.normCT)
      ax2.FontSize = fontsize * 0.85; 
      ax2.XTickLabelRotation = 0; 
      ax2.XTickLabel = {'0','','','0.5','','','1','','','1.5','','','2'};
      xlabel('normalized intensities with normalized SPM classes'); 
    else
      xlabel('CT intensities of SPM classes'); 
    end
    ax2.YTick = []; 
    if seg8t.isCTseg && ~opt.normCT 
      for ci = 1:3
        pl = plot([tismri.Tth(ci) tismri.Tth(ci)],[0 max(ylim)]); 
        pl.Color     = [ clscol(min(6,max(seg8t.lkp(ci))),:)];
      end
      for ci = 1:numel(tismri.iBonemn)
        pl = plot( repmat(tismri.iBonemn(ci),1,2) ,[0 max(ylim)]); 
        pl.Color = [ clscol(4,:)];
      end
      pl = plot( repmat( tis.report.fat ,1,2) ,[0 max(ylim)]); 
      pl.Color = [ clscol(5,:)];
      pl = plot( repmat( tis.report.muscle ,1,2) ,[0 max(ylim)]); 
      pl.Color = [ clscol(5,:)];
    else
      for ci = 1:numel(seg8t.mn)
        pl = plot([seg8t.mn(ci) seg8t.mn(ci)] / tis.WMth,[0 max(ylim)]); 
      end
      pl.Color     = [ clscol(min(6,max(seg8t.lkp(ci))),:) min(1,max(0,seg8t.mg(ci) * .7*sum(seg8t.lkp(ci)==seg8t.lkp).^.5))];
      pl.LineWidth = seg8t.mg(ci); %.^2 * 2; 
    end
    % backup?  plot([1 1],[0 max(ylim)],'-','Color',[.7 .7 .7]);
    if opt.method>0
      if numel(Yc)==6, lg = legend(flip({'GM','WM','CSF','bone','head'}),'box','off'); lg.Position(1) = .895; end 
    else
      legend(flip({'bone'}),'box','off'); 
    end
  end  
% ######################################
% horizontal boxplot for bones only with new measure as bars 
% ######################################



  %% == images ==
  if opt.report > 1
    spm_orthviews('Reset')
    pos = {[0.008 0.375 0.486 0.35]; [0.506 0.375 0.486 0.35]};
    % T1 + SPM segmentation
    if tis.weighting < 0 
      colormap( [[0 0.02 0.07]; repmat([0.05 0.15 .35],round(59/(crange+2)),1); repmat([ 0 .3 .6],round(59/(crange+2)),1); jet(59 - 2*round(59/(crange+2)))]);
      V0 = Vo; V0.dat(:,:,:) = single(.3*Yc{3} + .75*Yc{1} + 1*Yc{2} + 1.2*Yc{4} + 0.5*Yc{5} ); V0.dt(1) = 16;
      V0.pinfo  = repmat([1;0],1,size(Ym,3));
    elseif ~isempty(Ym)
    %colormap( [zeros(1,3); jet(59 - 2*round(59/(crange+2))); 0.5*ones(round(59/(crange+2)),3); ones(round(59/(crange+2)),3)])
      colormap( [[0 0.02 0.07]; repmat([0.05 0.15 .35],round(59/(crange+2)),1); repmat([ 0 .3 .6],round(59/(crange+2)),1); jet(59 - 2*round(59/(crange+2)))]);
      V0 = Vo; V0.dat(:,:,:) = single(Ym); V0.dt(1) = 16;
      V0.pinfo  = repmat([1;0],1,size(Ym,3));
    else
      colormap gray; 
      V0       = seg8t.image;
      Vc4      = spm_vol(spm_file(seg8t.image.fname,'prefix','c4'));
      Vc4.mat  = seg8t.Affine  * Vc4.mat; 
    end
    V0.mat    = seg8t.Affine * seg8t.image.mat; % Vo.mat;   \seg8t.tpm(1).mat
    hh0       = spm_orthviews('Image',spm_vol(fullfile(spm('dir'),'tpm','TPM.nii,1')),pos{1}); % avoid problems with high resolutions
    spm_orthviews('window',hh0,[0 10000]); % just make it black
    spm_orthviews('BB', [-85 -120 -90; 85 95 105]); % this has to be set in the low-resolution image
    hh0       = spm_orthviews('Image',V0,pos{1}); % add correct image after the other settings! 
    spm_orthviews('Reposition',[-25 0 0]);
    if tis.weighting < 0 
      spm_orthviews('Caption',hh0,sprintf('%s','CTseg'));
    else
      spm_orthviews('Caption',hh0,sprintf('%s',tis.weightingn));
    end
    if ~isempty(Ym)
      spm_orthviews('window',hh0,[0 1.2]); %
    else
      spm_orthviews('addtruecolourimage',hh0,Vc4,[0 0 0; 0.5 0.5 0.5; 1 0 0],0.4,1,0.1); 
      spm_orthviews('redraw');

      orthviewlegend = get(findobj(get(get(st.vols{1}.ax{1}.ax,'parent'),'children'),'Type','Image','Tag',''),'parent');
      orthviewlegend.YTick      = [0 .5 1];
      orthviewlegend.YTickLabel = {'BG','~Bone','Bone'};
      orthviewlegend.FontSize   =  orthviewlegend.FontSize * 0.85;
      
      try % does not work in headless mode without java
        figfs10 = [ findobj(fg,'FontSize',10); findobj(fg,'FontSize',9); findobj(fg,'FontSize',11); ...
          findobj(fg,'FontSize',8); findobj(fg,'FontSize',fontsize*.6);  findobj(fg,'FontSize',fontsize*.4); ]; 
        for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize * .8; end; end %#ok<TRYNC> 
        saveas(fg,Po.report,'png')
        for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize / .8; end; end %#ok<TRYNC> 
      catch
        cprintf('err','Cannot print report.')
      end
      
      return
    end
     %% 
    %  spm_orthviews('BB', [-90 -120 -80; 90 90 105]); %[[-78 -112 -70];[78 76 85]];
    %spm_orthviews('BB', [-90 -120 -80; 90 90 100] ./ [vx_vol;vx_vol]); % ######################## this is not corret but works ?! ... 
  
   
    
    % print cls overlay?
    if tis.weighting >= 0
      Vp0        = Vo; 
      Vp0.dt(1)  = 16;
      Vp0.mat    = seg8t.Affine * Vo.mat; 
      Vp0.pinfo  = repmat([1;0],1,size(Ym,3));
      switch 2
        case 1
          Vp0.dat(:,:,:) = single( (Yc{4}>.5) + (2*(Yc{5}>.5) + 0.5*(Yc{5}>.5)) ); 
          spm_orthviews('addtruecolourimage',hh0,Vp0,[0 0 0; 1 0 0; 1 1 0],0.4,2,0.1); 
        case 2
          Vp0.dat(:,:,:) = single( (Yc{4}>.5) + (Yc{4}>.5 & Ym>.7 ) + (3*(Yc{5}>.5)) + (Yc{5}>.5 & (Ym>.9)) ); 
          spm_orthviews('addtruecolourimage',hh0,Vp0,[0 0 0; 0 0.5 1; .5 0.8 1; .9 .5 0.5; 1 1 0 ],0.5,4,0); 
        case 5
          Vp0.dat(:,:,:) = single(2*(Yc{1}>.5) + 3*(Yc{2}>.5) + (Yc{3}>.5) + 4*(Yc{4}>.5) + 5*(Yc{5}>.5)); %.*(Ybone>.5) + 4*Yc{4}.*(Ybone>.5);
          spm_orthviews('addtruecolourimage',hh0,Vp0,[0 0 0; 0 0 1; 0 .8 0; 1 1 1; 1 0 0; 1 1 0],0.4,5,0); 
      end
    end
  
    %% print bone marrow
    if tis.weighting < 0
      V1 = Vo; V1.dat(:,:,:) = single((1000 + Ym)/2300); V1.dt(1) = 16;
    else
      Ybonemarrow2 = Ybonemarrow; Ybonemarrow2(isnan(Ybonemarrow2)) = 0; 
      V1 = Vo; V1.dat(:,:,:) = min(crange,Ybonemarrow2) + 2*(Ybonemarrow2>0) + single(Yc{1}>.5) + single(Yc{5}>.5) + 2*single(Yc{2}>0.5);
    end
    V1.pinfo  = repmat([1;0],1,size(Ym,3));
    V1.mat    = seg8t.Affine * V1.mat; 
    V1.dt(1)  = 16;
    hh1       = spm_orthviews('Image',V1,pos{2}); 
    spm_orthviews('Interp',1);
    if tis.weighting < 0
      spm_orthviews('window',hh1,[0 1.2]);
      spm_orthviews('Caption',hh1,sprintf('%s',tis.weightingn));
    else
      spm_orthviews('window',hh1,[0 crange+2]);
      spm_orthviews('Caption',hh1,'bonemarrow + GM/HD + 2*WM');
    end
    spm_orthviews('Reposition',[-25 0 0]);
    spm_orthviews('redraw'); 
   
    
    % this is replaced by spm_orthviews('redraw'); 
    if exist('orthviewlegend','var')
      if opt.method > 0
        orthviewlegend = get(findobj(get(get(st.vols{1}.ax{1}.ax,'parent'),'children'),'Type','Image','Tag',''),'parent');
        orthviewlegend.YTickLabel = {'BG','lBone','hBone','muscle','fat'};
        orthviewlegend.FontSize   =  orthviewlegend.FontSize * 0.85;
      else
        orthviewlegend = get(findobj(get(get(st.vols{1}.ax{1}.ax,'parent'),'children'),'Type','Image','Tag',''),'parent');
        orthviewlegend.YTick      = {0 1};
        orthviewlegend.YTickLabel = {'BG','Bone'};
        orthviewlegend.FontSize   =  orthviewlegend.FontSize * 0.85;
      end
    end

    %% bone surface 
    if isfield(Po,'central') && ~isempty(St)
      for idi = 2:3 %1:2
        spm_orthviews('AddContext',idi);  
        spm_ov_mesh('display',idi,Po.central );
        % apply affine transformation
        V = (seg8t.Affine * ([st.vols{idi}.mesh.meshes(end).vertices,...
             ones(size(st.vols{idi}.mesh.meshes(end).vertices,1),1)])' )';
        V(:,4) = [];
        st.vols{idi}.mesh.meshes = subsasgn(st.vols{idi}.mesh.meshes,struct('subs','vertices','type','.'),single(V));
        % change line style
        hM = findobj(st.vols{idi}.ax{1}.cm,'Label','Mesh');
        UD = get(hM,'UserData');  
        UD.width = 1;
        if tis.weighting == -1
          UD.style = {'w-'};  
        else
          if idi==2, UD.style = {'r-'}; else, if tis.boneIntType>1, UD.style = {'w-'}; else, UD.style = {'r-'}; end; end
          %UD.style = {'r-'};  
        end
        set(hM,'UserData',UD); clear UD hM
        warning('off','MATLAB:subscripting:noSubscriptsSpecified');
        spm_ov_mesh('redraw',idi);
      end
    end
  end
 




  %% == surface ==
  clear h; 
  if ~isempty(Stm) && opt.report > 1
    % thickness map
    hCS{1} = subplot('Position',[0.015 0.20 0.23 0.15],'Parent',fg,'visible','off');  sview{1} = 'l';
    hCS{2} = subplot('Position',[0.255 0.20 0.23 0.15],'Parent',fg,'visible','off');  sview{2} = 'r';
    hCS{3} = subplot('Position',[0.015 0.01 0.23 0.19],'Parent',fg,'visible','off');  sview{3} = 't';
    hCS{4} = subplot('Position',[0.255 0.06 0.23 0.14],'Parent',fg,'visible','off');  sview{4} = 'p';
    hCS{5} = subplot('Position',[0.255 0.03 0.23 0.02],'Parent',fg,'visible','off'); 
    
    imat = spm_imatrix(seg8t.Affine); Rigid = spm_matrix([imat(1:6) ones(1,3)*mean(imat(7:9)) 0 0 0]); clear imat;
    V = (Rigid * ([Stm.vertices, ones(size(Stm.vertices,1),1)])' )'; V(:,4) = []; Stm.vertices = V;

    for ci = 1:4
      h{ci} = cat_surf_render2(Stm,'parent',hCS{ci}); 
      cat_surf_render2('Clim',h{ci},[0 20]); 
      switch lower(sview{ci})
        case {'r'},  cat_surf_render('view',h{ci},[  90   0]); 
        case {'l'},  cat_surf_render('view',h{ci},[ -90   0]);  
        case {'t'},  cat_surf_render('view',h{ci},[   0  90]); 
        case {'p'},  cat_surf_render('view',h{ci},[   0   0]); 
      end
    end
    
    cb = cat_surf_render2('Colorbar',h{4});
    cb.colourbar.Location = 'South';
    cb.colourbar.Position = [0.275 0.015 0.2 0.005];
    text(0.5,0.015,'Bone thickness (tBone)','fontsize',fontsize-1,'FontName',fontname,...
        'HorizontalAlignment','center','Parent',hCS{5}); 

  end

  if ~isempty(Si) && opt.report > 1
    hCS{1} = subplot('Position',[0.5 + 0.015 0.20 0.23 0.15],'Parent',fg,'visible','off');  sview{1} = 'l';
    hCS{2} = subplot('Position',[0.5 + 0.255 0.20 0.23 0.15],'Parent',fg,'visible','off');  sview{2} = 'r';
    hCS{3} = subplot('Position',[0.5 + 0.015 0.01 0.23 0.19],'Parent',fg,'visible','off');  sview{3} = 't';
    hCS{4} = subplot('Position',[0.5 + 0.255 0.06 0.23 0.14],'Parent',fg,'visible','off');  sview{4} = 'p';
    hCS{5} = subplot('Position',[0.5 + 0.255 0.03 0.23 0.02],'Parent',fg,'visible','off'); 

    imat = spm_imatrix(seg8t.Affine); Rigid = spm_matrix([imat(1:6) ones(1,3)*mean(imat(7:9)) 0 0 0]); clear imat;
    V = (Rigid * ([Si.vertices, ones(size(Si.vertices,1),1)])' )'; V(:,4) = []; Si.vertices = V;

    for ci = 1:4
      h{ci} = cat_surf_render2(Si,'parent',hCS{ci}); 
      if tis.weighting < 0 
        cat_surf_render2('Clim',h{ci},[500 1500]); 
      else
        cat_surf_render2('Clim',h{ci},[0 crange]); 
      end
      switch lower(sview{ci})
        case {'r'},  cat_surf_render('view',h{ci},[  90   0]); 
        case {'l'},  cat_surf_render('view',h{ci},[ -90   0]);  
        case {'t'},  cat_surf_render('view',h{ci},[   0  90]); 
        case {'p'},  cat_surf_render('view',h{ci},[   0   0]); 
      end
    end

    if 1
      cb = cat_surf_render2('Colorbar',h{4});
      cb.colourbar.Location = 'South';
      cb.colourbar.Position = [0.5 + 0.275 0.015 0.2 0.005];
      text(0.5,0.015,'Bone intensity (iBone)','fontsize',fontsize-1,'FontName',fontname,...
        'HorizontalAlignment','center','Parent',hCS{5}); 

    else
      %% colormap
      if strcmpi(spm_check_version,'octave')
        axes('Position',[0.965 0.03 0.01 0.28],'Parent',fg); image(flip(121:1:120+surfcolors)','Parent',cc{4});
        cc{4} = gca; 
      else
        cc{4} = axes('Position',[0.965 0.03 0.01 0.28],'Parent',fg); image(flipud(121:1:120+surfcolors)','Parent',cc{4});
      end
      
      %% histogram line
      if strcmpi(spm_check_version,'octave')
        axes('Position',[0.936 0.03 0.03 0.28],'Parent',fg,'Visible', 'off','tag', 'cat_surf_results_hist', ...
        'xcolor',fontcolor,'ycolor',fontcolor);
        cc{5} = gca; 
      else
        cc{5} = axes('Position',[0.936 0.03 0.03 0.28],'Parent',fg,'Visible', 'off','tag', 'cat_surf_results_hist', ...
          'xcolor',fontcolor,'ycolor',fontcolor);
      end
      side  = hSD{1}.cdata;
      [d,h] = hist( side(~isinf(side(:)) & ~isnan(side(:)) &  side(:)<6 & side(:)>0) ,  hrange); %#ok<HIST> 
      d = d./numel(side);
      d = d./max(d);
      text(cc{5},'test')
      
      % print histogram
      hold(cc{5},'on');  
      for bi = 1:numel(d)
        b(bi) = barh(cc{5},h(bi),-d(bi),boxwidth); 
        set(b(bi),'Facecolor',cmap3(bi,:),'Edgecolor',fontcolor); 
      end
      ylim([0,20]); xlim([-1 0]);
    end
  end


  %
  
  colormap( [[0 0.02 0.07]; repmat([0.05 0.15 .35],round(59/(crange+2)),1); repmat([ 0 .3 .6],round(59/(crange+2)),1); jet(59 - 2*round(59/(crange+2)))]);
    

  try % does not work in headless mode without java
    figfs10 = [ findobj(fg,'FontSize',fontsize+1); findobj(fg,'FontSize',fontsize); findobj(fg,'FontSize',fontsize-1);  ...
      findobj(fg,'FontSize',fontsize*0.85); findobj(fg,'FontSize',fontsize*.6);  findobj(fg,'FontSize',fontsize*.4); ]; 
    for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize * .75; end; end %#ok<TRYNC> 
    saveas(fg,Po.report,'png')
    for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize / .75; end; end %#ok<TRYNC> 
  end

end
%=======================================================================
function [Theader,Tline,Tavg, Cheader, MAfn, matm, mmatm]   = prepare_print(P,opt)
%prepare_print. Create table elements for the in-line report.
  if 0
    MAfn    = {'Tw','Thbg','Tfat','Tbone','Tres','Tcnr',      'BG',  'CSF','GM',   'muscle','fat',   'bone','marrow',   'MBI' , 'BT' , 'HDT' }; % 'Pll', Tw, Tfat, Thg, Tcnr
    MAffn   = {'s' ,'s'   ,'s'   ,'s'    ,'f'   ,'f'   ,      'f' ,  'f'  ,'f' ,   'f'     ,'f'  ,   'f'   ,'f'     ,   'f'   , 'f'  , 'f'   }; % 'e'  ,
    MAfsep  = [0   ,0     ,0     ,0      ,0     ,0     ,      1   ,  0    ,0   ,   1       ,0    ,   1     ,0       ,   0     , 0    , 0     ]; % 0    ,
  else
    MAfn    = {'Tw','Thbg','Tfat','Tbone','Tres','Tcnr',      'BG',  'CSF','GM',   'Tskull','Thead',  'MED' , 'MBI'  , 'BT' , 'HDT' }; % 'Pll', Tw, Tfat, Thg, Tcnr
    MAffn   = {'s' ,'s'   ,'s'   ,'s'    ,'f'   ,'f'   ,      'f' ,  'f'  ,'f' ,   'f'     ,'f'    ,  'f'   , 'f'    , 'f'  , 'f'   }; % 'e'  ,
    MAfsep  = [0   ,0     ,0     ,0      ,0     ,0     ,      1   ,  0    ,0   ,   1       ,0      ,  1     ,  0     , 0    , 0 ]; % 0    ,
  end
  
  Cheader = {'scan'};
  Theader = sprintf(sprintf('%%%ds:',opt.snspace(1)-1),'scan');
  Tline   = sprintf('%%5d) %%%ds:',opt.snspace(1)-8);
  Tavg    = sprintf('%%%ds:',opt.snspace(1)-1);
  for fi = 1:numel(MAfn)
    Cheader = [Cheader MAfn{fi}]; 
    AMfni   = strrep( MAfn{fi} ,'_','');
    if MAfsep(fi), sep = ' |'; else, sep = ''; end
    Theader = sprintf(sprintf('%%s%%s%%%ds' ,opt.snspace(2)),Theader, sep, AMfni ); 
    switch MAffn{fi}
      case {'s','d'}
        Tline   = sprintf('%s%s%%%d%s', Tline , sep, opt.snspace(2), MAffn{fi});
        Tavg    = sprintf('%s%s%%%d%s', Tavg  , sep, opt.snspace(2), MAffn{fi});
      otherwise
        Tline   = sprintf('%s%s%%%d.%d%s',Tline , sep, opt.snspace(2), opt.snspace(3) .* (MAffn{fi}=='f'),MAffn{fi});
        Tavg    = sprintf('%s%s%%%d.%d%s',Tavg  , sep, opt.snspace(2), opt.snspace(3) .* (MAffn{fi}=='f'),MAffn{fi});
    end
  end
  % add time 
  Tline   = sprintf('%s%%4.0fs',Tline);
  
  % main result table
  matm    = num2cell(nan(numel(P),numel(MAfn)));
  mmatm   = 10.5*ones(numel(P),numel(MAfn));
  
  % print title
  if opt.verb
    fprintf('\nBone Preprocessing:\n');
    if opt.verb 
      methodstr = {'SPMmat8','MedBoneSurf','rMedBoneSurf'};
      reportstr = {'Table','Table + Volumes','Table + Volumes + Surfaces'};
      fprintf('  Method:   %d (%s)\n', opt.method, methodstr{opt.method+1});
      fprintf('  Report:   %d (%s)\n', opt.report, reportstr{opt.report});
      if opt.method>0
        fprintf('  Reduce:   %d (%s)\n', opt.reduce, reportstr{opt.reduce});
      end
      % parameter ######################
      %  - method id + name it 
      %  - report option + name it (only tables, imags, full)
      %  - warnings/information for atypical settings
      % ################################
    end
    fprintf('\n%s\n%s\n',  Theader,repmat('-',size(Theader)));  
  end
end
%=======================================================================
function [ seg8t, tis, matm, mmatm, vx_vol ]          = get_segmat(out,MAfn)   % TODO: (1) lkp tests; (2) skull-stripping/defacing tests; (3) TIV estimate  
%get_segmat. Get and evaluate SPM preprocessing structure.
%
%  [seg8t, tis, tis, matm, vx_vol] = get_spm8(P,MAfn)
%
%   P       .. one image
%   MAfn    .. fieldnames for matm
%
%   seg8t   .. spm8 structure without larger fields
%   tis     .. measures based on SPM tissue thresholds and 
%              basic information based on image and tissue properties
%   matm    .. line-wise output (print value)
%   mmatm   .. line-wise output (mark  value)
% 
% Comments: 
% * SPM tissue peaks for [GM WM CSF bone head background]
%   - default is [1 1 2 3 4 2]: 
%     . in T1 the lower CSF peak is the right one whereas the other one is PVE to GM or in the best case meninges?  
%       > would need a seperate class?
%     . the lower CSF threshold trend to overestimation in younger subjects - PVE?
%       > the CSF value is less robust for normalization 
%     . in mixed tissues, higher mg result often in higher vr that bias the mn  


% load SPM mat 
  seg8t          = load(out.P.seg8); 
  if out.CTseg
    seg8t.dat.model.gmm = rmfield(seg8t.dat.model.gmm,{'T','Sig'});
    % #### this is not fully working and the classe values are strange ...
    seg8t.lkp     = seg8t.sett.gmm.mg_ix;
    seg8t.mn      = seg8t.dat.model.gmm.m; 
    seg8t.mg      = seg8t.dat.model.gmm.gam'; 
    seg8t.vr      = seg8t.dat.model.gmm.W; 
     
    tmp           = spm_load_priors8(ps_fullfile(spm('dir'),'TPM','TPM.nii'));
    seg8t.tpmA    = tmp;
    seg8t.tpm     = rmfield(tmp.V,'private');
    seg8t.Affine  = eye(4);
    seg8t.image   = spm_vol(out.P.org);
    seg8t.isCTseg = 1; 
    seg8t         = rmfield(seg8t,'dat');
    seg8t.sett    = rmfield(seg8t.sett,'B');
  else
    seg8t         = rmfield(seg8t,{'Twarp','Tbias','MT'}); % store this one for later 
    seg8t.isCTseg = 0; 
    
    if numel( seg8t.image ) > 1
      iid         = contains(out.P.org,{seg8t.image.fname}); 
      seg8t.image = seg8t.image(iid);
      seg8t.mn    = seg8t.mn(iid,:);
      seg8t.vr    = seg8t.vr(iid,iid,:);
    end
  end
  vx_vol         = sqrt(sum(seg8t.image(1).mat(1:3,1:3).^2));
  


  % check for problems and skip in worst case
  if numel(seg8t.tpm) ~= 6 
    cat_io_cprintf('err','ERROR: Only 6 class models are supported yet!\n');
    return
  end
  % #####################
  % check number of Gaussian peak per class ? 
  
  

 
  % == evaluate SPM mat entries to name some basic images properties ==
  % SPM main tissue thresholds
  tis.help.seg8o      = 'SPM seg8 main tissue intensity (mn of max mg in lkp)';
  tis.help.seg8ov     = 'SPM seg8 main tissue variance (vr of max mg in lkp)';
  tis.help.seg8n      = 'SPM main tissue intensity (mn of max mg in lkp) normalized by the WM'; % normalized by WM
  tis.help.seg8nv     = 'SPM main tissue variance (vr of max mg in lkp) normalized by the WM';
  tis.help.seg8con    = 'mimimum brain tissue contrast in SPM seg8t'; 
  tis.help.seg8conn   = 'mimimum brain tissue contrast in SPM seg8t normalized by the WM'; 
  tis.help.seg8CNR    = 'Noise estimate as minimum of the WM and CSF variance in percent (similar to BWP).';
  tis.help.WMth       = 'SPM WM tisse intensity.'; 
  tis.help.res_vx_vol = 'Image voxel resolution in mm.'; 
  tis.help.res_RES    = 'RMS voxel resolution.'; 
  tis.seg8o           = nan(1,6);
  tis.seg8ov          = nan(1,6);
  for ci = 1:max(seg8t.lkp) 
    % The SPM Gaussian's seams to be unsortet and sorting based on the mean
    % value or the variance would be useful 
    sortvar = 'vr';
    [~,sorti] = sort( seg8t.(sortvar)(seg8t.lkp==ci) ); 
    var = seg8t.mn( seg8t.lkp==ci ); tis.seg8mns( seg8t.lkp==ci ) = var(sorti); 
    var = seg8t.mg( seg8t.lkp==ci ); tis.seg8mgs( seg8t.lkp==ci ) = var(sorti); 
    var = seg8t.vr( seg8t.lkp==ci ); tis.seg8vrs( seg8t.lkp==ci ) = var(sorti); 
    
    % How to average values ... well, when we are interested in changes of
    % intensities in bone(marrow) and the propostion of fat then weighted  
    % averaging should be fine (will depend on subject and protocol).
    if 0 % use major class, ie. high proposion and low variance 
      tis.seg8o(ci)   = seg8t.mn( (seg8t.lkp==ci) & seg8t.mg' == max( seg8t.mg' .* (seg8t.lkp==ci)) );
      tis.seg8ov(ci)  = seg8t.vr( (seg8t.lkp==ci) & seg8t.mg' == max( seg8t.mg' .* (seg8t.lkp==ci)) );
    else % mix peak as average class intensity ... 
      tis.seg8o(ci)   = mean(seg8t.mg(seg8t.lkp==ci)' .* seg8t.mn(seg8t.lkp==ci),2);
      tis.seg8ov(ci)  = mean(seg8t.mg(seg8t.lkp==ci)' .* shiftdim(seg8t.vr(seg8t.lkp==ci),4),2);
    end
    if 0 %ci == 3 && sum(seg8t.lkp==ci)>1 % CSF subcase to avoid PVE in T1 images
      if mean( seg8t.mn(seg8t.lkp==ci) )  <  mean( seg8t.mn(seg8t.lkp<=2) ) % low CSF: CSF < WM+GM
        csfmn = seg8t.mn( seg8t.lkp==ci ); csfvr = seg8t.vr( seg8t.lkp==ci );
        [tis.seg8o(ci),csfmnmin] = min( csfmn );
        tis.seg8ov(ci) = csfvr(csfmnmin);
  
        tis.help.report_meninges = 'In case of multiple CSF Gaussians, the larger one in T1 often presents the intensity of meninges.'; 
        tis.meninges = mean( csfmn( setdiff(1:numel(csfmn),csfmnmin) ));
      else
        tis.help.report_meninges = 'In case of only one CSF Gaussian, the value of meninges should be in the middle between CSF and GM.'; 
        tis.meninges = mean( [tis.seg8ov(1) tis.seg8ov(3)] ); 
        % do not know enough about T2/PD case
      end
    end
  end

  % test if an image is a CT based on the intensities
  % low tissue contrast, but high bone, and low background
  %isCT = out.CTseg; % ###############  not working
       %  (tis.seg8o(1:3) > 0   & tis.seg8o(1:3) < 50) & ...  
       %  (tis.seg8o(4) > 600   & tis.seg8o(4) < 2000) & ...
       %  (tis.seg8o(5) > -50   & tis.seg8o(5) < 50)   & ...
       %  (tis.seg8o(6) > -2000 & tis.seg8o(6) < -600 );
 
  if out.CTseg
    % test CT condition 
    tis.WMth          = 40; 
    tis.seg8o         = [40 30 20 1024 0 -1024];
    tis.seg8n         = [40 30 20 1024 0 -1024] / 2000 +1;
    tis.seg8nv        = [40 30 20 1024 0 -1024];
    tis.seg8con       = nan; 
    tis.seg8conr      = nan; 
    tis.seg8CNR       = nan; 
  else
    tis.WMth          = tis.seg8o(2);
    tis.seg8n         = tis.seg8o  ./ tis.WMth; % normalized by WM
    tis.seg8nv        = tis.seg8ov ./ tis.WMth;
    tis.seg8con       = min(abs(diff(tis.seg8o)));
    tis.seg8conr      = min(abs(diff(tis.seg8n)));
    tis.seg8CNR       = min( seg8t.vr(:) ./ tis.WMth ./ (seg8t.lkp(:)==2 | seg8t.lkp(:)==3  | seg8t.lkp(:)==6) ) .* tis.seg8conr * 3; 
  end
  tis.res_vx_vol    = vx_vol; 
  tis.res_RES       = mean(vx_vol.^2).^.5;
  
  tis.spm_bone_con  = sum(seg8t.mn( seg8t.lkp == 3 )' .* ...
                      ( seg8t.mg( seg8t.lkp == 3 ) == max(seg8t.mg( seg8t.lkp == 3 ) ) ));
            
  if max(seg8t.lkp) > 6
      tis.spm_bone_med1 = seg8t.mn( seg8t.lkp == 4 ) / tis.spm_bone_con;
      tis.spm_head_med1 = seg8t.mn( seg8t.lkp == 5 ) / tis.spm_bone_con;
  else
      % default case for 6 TPM classes
      tis.spm_con = sum(seg8t.mn( seg8t.lkp == 3 )' .* ...
        ( seg8t.mg( seg8t.lkp == 3 ) == max(seg8t.mg( seg8t.lkp == 3 ) ) ));
      [bone,bonei] = sort(seg8t.mn( seg8t.lkp == 4 ) / tis.spm_con); % order intensities within class 4 (normalized by CSF)
      bonevr = seg8t.vr( seg8t.lkp == 4 ) / tis.spm_con;
      if sum(seg8t.lkp == 4) == 3
        tis.spm_bone_med1   = bone(3); % highest value
        tis.spm_bone_med2   = bone(1); 
        tis.spm_bone_med3   = bone(2); 
        tis.spm_bone_medvr1 = bonevr(bonei(3)); % highest value
        tis.spm_bone_medvr2 = bonevr(bonei(1)); 
        tis.spm_bone_medvr3 = bonevr(bonei(2)); 
      elseif sum(seg8t.lkp == 4) == 2 % CTseg
        tis.spm_bone_med1 = bone(1); % highest value
        tis.spm_bone_med2 = bone(2); 
        tis.spm_bone_medvr1 = bonevr(bonei(1)); % highest value
        tis.spm_bone_medvr2 = bonevr(bonei(2)); 
      else
        tis.spm_bone_med1   = bone(1); % highest value
        tis.spm_bone_med2   = bone(1); 
        tis.spm_bone_medvr1 = bonevr(bonei(1)); % highest value
        tis.spm_bone_medvr2 = bonevr(bonei(1)); 
      end
  end


  % image weighting 
  tis.help.weighting = 'Image weigthing based on SPM seg8t intensities (0=PDw; 1=T1w; 2=T2w; 3=MTw; 4=IRw, -1=CT).'; 
  if out.CTseg
      tis.weighting  = -1;
      tis.weightingn = 'CT';
  else
    if     tis.seg8n(3) < tis.seg8n(1)  &&  tis.seg8n(1) < tis.seg8n(2)      % T1: CSF < GM < WM 
      tis.weighting  = 1; 
      tis.weightingn = 'T1w';
    elseif tis.seg8n(3) > tis.seg8n(1)  &&  tis.seg8n(1) > tis.seg8n(2)      % T2: CSF > GM > WM
      tis.weighting  = 2;
      tis.weightingn = 'T2w';
    elseif tis.seg8o(3) < 0  && is.seg8o(2) < 2                              % MT: negative CSF values and not to high WM values
      tis.weighting  = 3;
      tis.weightingn = 'MTw';
    elseif tis.seg8o(6) > tis.seg8o(2) && tis.seg8nv(6) < .3                 % high int - low var
      tis.weighting  = 4;
      tis.weightingn = 'IRw'; % inverse recovery
    else
      tis.weighting  = 0; 
      tis.weightingn = 'PDw';
    end
  end


  % background types 
  tis.help.highBG = ['Intensity of the background based on the SPM seg8t intensities ' ...
    '(0=low,classical MRI; 1=high int low var, eg. MP2Rage; 2=high int low var, eg. IR; 3=mid int high var, eg. MT)'];
  if     tis.seg8o(6) > tis.seg8o(2) && tis.seg8nv(6) > .3 % high int - high var
    tis.highBG     = 1;
    tis.highBGn    = 'high'; % intensity, high variance background (eg. uncorrected MP2Rage)';
  elseif tis.seg8o(6) > tis.seg8o(2) && tis.seg8nv(6) < .3 % high int - low var
    tis.highBG     = 2;
    tis.highBGn    = 'high2'; % intensity, low variance background (eg. inverse recovery)';
  elseif tis.seg8o(6) < tis.seg8o(3) 
    tis.highBG     = 0;
    tis.highBGn    = 'low'; % intensity, low variance background (eg. classical MRI)';
  else
    tis.highBG     = 3;
    tis.highBGn    = 'mid'; % intensity, high variance background (eg. MT)';
    tis.weighting  = 3;
    tis.weightingn = 'MTw';
  end


  % fat suppression
  tis.help.minBone     = 'Minimum Gaussian of the bone tissue class with more as 5% normalized by the WM intensity. ';
  tis.help.headFatType = 'Protocol intensity type of the head based on SPM seg8t values (0-low[<CSV], 1-mid[<WM], 2-[>WM]). ';
  tis.help.boneIntType = 'Protocol intensity type of the bone based on SPM seg8t values (0-low[<CSV], 1-mid[<WM], 2-[>WM]). ';
  tis.minHead = min(    seg8t.mn (seg8t.lkp==5 & seg8t.mg'>0.05)) ./ tis.seg8o(2); 
  tis.medHead = median( seg8t.mn (seg8t.lkp==5 & seg8t.mg'>0.05)) ./ tis.seg8o(2); 
  tis.maxHead = max(    seg8t.mn (seg8t.lkp==5 & seg8t.mg'>0.05)) ./ tis.seg8o(2); 
  tis.minBone = min(    seg8t.mn (seg8t.lkp==4 & seg8t.mg'>0.05)) ./ tis.seg8o(2);
  tis.medBone = median( seg8t.mn (seg8t.lkp==4 & seg8t.mg'>0.05)) ./ tis.seg8o(2);
  tis.maxBone = max(    seg8t.mn (seg8t.lkp==4 & seg8t.mg'>0.05)) ./ tis.seg8o(2);
  if out.CTseg
    tis.headFatType   = 0; 
    tis.headFatTypen  = '-';
    tis.boneIntType   = 0; 
    tis.boneIntTypen  = '-'; 
    % bone and head values of CTseg are not realy useful :/ 
  else
    if tis.maxHead > 1.2
      tis.headFatType  = 2; 
      tis.headFatTypen = 'high'; 
    elseif tis.maxHead < tis.seg8n(3) 
      tis.headFatType  = 0; 
      tis.headFatTypen = 'low'; 
    else
      tis.headFatType  = 1; 
      tis.headFatTypen = 'mid'; 
    end
    if tis.headFatType > 1  &&  tis.maxBone > 1
      tis.boneIntType   = 2; 
      tis.boneIntTypen  = 'high'; 
    elseif tis.headFatType < 2  &&  tis.minBone < 0.1
      tis.boneIntType   = 0; 
      tis.boneIntTypen  = 'low'; 
    else
      tis.boneIntType   = 1; 
      tis.boneIntTypen  = 'mid'; 
    end
  end  
    




  % skull-stripped ? 
  % #######################################################
  % create error / continue


  % defacing ?
  % deearing ?
  % >> affected area ? 
  %    surface percentage betwee cutted (0-values) and uncutted / full 
  %    background area
  % #######################################################
    
  
  % TIV estimate ? 
  % ############################
  % evaluate Affine value? > maybe with some fix TPM TIV value?
  % evaluate log-likelihood?



  if 0
  %% just an internal overview of SPM paraemter
    WMth = seg8t.mn( (seg8t.lkp==2) &  seg8t.mg' == max( seg8t.mg' .* (seg8t.lkp==2) )); % main WM intensity 
    fprintf('Tissue table [class; propotion; WM normed intensity; WM normed var; normed var]:\n')
    disp( [seg8t.lkp; seg8t.mg';  
      seg8t.mn / WMth; shiftdim(seg8t.vr,1)  ./ WMth; % normalized by WM intensity 
      shiftdim(seg8t.vr,1) ./ WMth  .* (seg8t.mn / WMth);  ] ) 
     
  end




  % report tissue types
  % - help
  tis.help.report        = 'Different measures used in the command line report.';
  tis.help.report_Tw     = 'tis.weightingn';
  tis.help.report_Thbg   = 'tis.highBGn'; 
  tis.help.report_Tfat   = 'tis.headFatTypen'; 
  tis.help.report_Tbone  = 'tis.boneFatTypen'; 
  tis.help.report_Tres   = 'tis.res_RES'; 
  tis.help.report_Tcnr   = 'tis.seg8CNR';
  tis.help.report_Pll    = 'seg8t.ll'; 
  tis.help.report_BG     = 'tis.seg8n(6)';
  tis.help.report_CSF    = 'tis.seg8n(3)';
  tis.help.report_GM     = 'tis.seg8n(1)';
  tis.help.report_WM     = 'tis.seg8n(2)';
  tis.help.report_muscle = 'normalized protocol-selected values of tis.*Head.'; 
  tis.help.report_fat    = 'normalized protocol-selected values of tis.*Head.'; 
  tis.help.report_bone   = 'tis.minBone';  
  tis.help.report_marrow = 'tis.maxBone'; 
  tis.help.report_MBI    = 'bone mean intensity (tis.maxBone - tis.report.bone) / (tis.report.fat - tis.report.bone) * 2'; 
  
  % - measures
  tis.report.Tw     = tis.weightingn; 
  tis.report.Thbg   = tis.highBGn; 
  tis.report.Tskull = tis.seg8n(4); 
  tis.report.Thead  = tis.seg8n(5); 
  tis.report.Tfat   = tis.headFatTypen; 
  tis.report.Tbone  = tis.boneIntTypen; 
  tis.report.Tres   = tis.res_RES; 
  tis.report.Tcnr   = tis.seg8CNR;
  %tis.report.Pll    = seg8t.ll; 
  % brain tissue peaks 
  tis.report.BG     = tis.seg8n(6);
  tis.report.CSF    = tis.seg8n(3);
  tis.report.GM     = tis.seg8n(1);
  tis.report.WM     = tis.seg8n(2);
  % head tissues peaks 
  if tis.weighting == 1 % T1 
    musth = (seg8t.lkp==5) & (seg8t.mg'>0.05) & (seg8t.mn > tis.seg8o(3)) & (seg8t.mn < tis.seg8o(2)); % head mn value between CSF and WM
    fatth = (seg8t.lkp==5) & (seg8t.mg'>0.05) & (seg8t.mn > tis.seg8o(2)); 
    if sum(musth)>0,  tis.report.muscle = min(seg8t.mn(musth)) / tis.seg8o(2); else; tis.report.muscle = tis.medHead; end
    if sum(fatth)>0,  tis.report.fat    = max(seg8t.mn(fatth)) / tis.seg8o(2); else; tis.report.fat    = tis.maxHead; end
  else
     tis.report.muscle = tis.medHead; 
     tis.report.fat    = tis.minHead; 
  end
  % bone tissues peaks
  tis.report.bone   = tis.minBone;  
  tis.report.marrow = tis.maxBone; 
  % ###################### here we need to account the percentage values ;-)  
  if tis.boneIntType == 1
    tis.report.MBI    = (tis.maxBone - tis.report.bone) / (tis.report.fat - tis.report.bone) * 2; 
  else
    tis.report.MBI    = 1 - ( (tis.maxBone - tis.report.bone) / (tis.report.fat - tis.report.bone) * 2); 
  end
  tis.report.MED      = nan;  
  tis.report.BT       = nan; 
  tis.report.HDT      = nan; 
  
  % create one line for the output table depending on the MAfn defintion
  for fni = 1:numel(MAfn)
    matm{1,fni}  = tis.report.(MAfn{fni});
    if isnumeric(tis.report.(MAfn{fni}))
      mmatm(1,fni) = tis.report.(MAfn{fni});
    else
      mmatm(1,fni) = 1; 
    end
  end

end
%=======================================================================
function [Vo, Yo, Yc, Ya, Ymsk, Ym, Affine, RES, BB]        = loadMRI(P,opt,seg8t,tis,cls,bd) % TODO: (1) tissue-base intensity scaling; (2) use/eval affreg?
% loadMRI. Load MRI and segmentation.
%
%  [Vo,Yo,Yc,Ya,Ym, Affine] = loadMRI(P,Pa)
%
%  Po,Vo,Yo .. original image
%  Pc,Vc,Yo .. segment class images (cell)
%  Pa,Va,Ya .. atlas image
%

%opt.affreg, 1:5, opt.reslim, 25 ,opt); 
  
  if ~exist('cls'   ,'var'), cls    = 1:5; end 
  if ~exist('bd'    ,'var'), bd     = 25;  end 

  % get bias corrected original image
  Vo = spm_vol(P.bc);
  Yo = single(spm_read_vols(Vo));
  
  % load segmentation 
  Pc = cell(1,6); Yc = cell(1,6); Yc{6} = ones(Vo.dim,'single');
  for ci = cls % 1:5
    if seg8t.isCTseg % CTseg 
      Pc{ci}  = fullfile(P.orgpp,sprintf('c%02d%s%s',ci,P.ppff(4:end),P.ee));
    else
      Pc{ci}  = fullfile(P.orgpp,sprintf('c%d%s%s',ci,P.orgff,P.ee));
    end
    Vc(ci)  = spm_vol(Pc{ci}); 
    Yc{ci}  = single(spm_read_vols(Vc(ci)));
    Yc{6}   = Yc{6} - Yc{ci}; 
  end

  % create a linear intensity normalized image
  % .. unclear side effects ... and the histogram is not looking nice ?
  if 0 %~isempty(Pa)
    minimg = min( Yo(:) ); 
    maximg = max( Yo(:) );
    mintis = min( tis.seg8o );
    maxtis = max( tis.seg8o(2)*1.5 , tis.maxHead*tis.seg8o(2) ); 
    switch tis.weighting 
      case 1 % T1
        isc   = 1; 
        T3th  = [minimg mintis tis.seg8o(3) tis.seg8o(1) tis.seg8o(2) tis.seg8o(2)+diff(tis.seg8o(1:2)) maxtis maximg]; 
        T3thx = [0 0.05 1 2 3 4 5 6]; 
        T3th  = interp1(T3th ,1:1/isc:numel(T3th )*isc,'pchip'); T3th  = smooth(T3th ,16*isc); %spm_smooth(T3th ,T3th ,.2*isc)
        T3thx = interp1(T3thx,1:1/isc:numel(T3thx)*isc,'pchip'); T3thx = smooth(T3thx,16*isc); %spm_smooth(T3thx,T3thx,.2*isc)
        Ym    = cat_main_gintnorm(Yo,struct('T3th',T3th,'T3thx',T3thx));
      case 2 % T2
        Ym = cat_main_gintnorm(Yo,struct('T3th',[minimg mintis tis.seg8o(2) ...
          tis.seg8o(1) tis.seg8o(3) maxtis maximg],'T3thx',[0 0.05 1 2 3 5 6]));
      case 3
        Ym = (Yo - min([ max([-.5 minimg]) tis.seg8o ])) / (tis.seg8o(2) - min([ max([-.5 minimg]) tis.seg8o ]));
      otherwise
        Ym = (Yo - min([ 0 tis.seg8o(3),tis.seg8o(end)])) / (tis.seg8o(2) - min([ 0 tis.seg8o(3),tis.seg8o(end)]));
    end
  else
    if tis.weighting == 2 % MT
      Ym = (Yo - min([-.5 tis.seg8o ])) / (tis.seg8o(2) - min([-.5 tis.seg8o ]));
    elseif tis.weighting == -1 % CT
      if opt.normCT 
        Ym = (Yo - min( tis.seg8o )) / max(tis.seg8o(:) - min(tis.seg8o(:)));
      else
        Ym = Yo; 
      end
    else
      Ym = (Yo - min( tis.seg8o )) / (tis.seg8o(2) - min([tis.seg8o(3),tis.seg8o(end)]));
    end
  end
  


  % == do affine registration ==
  % ##################### & opt.refine ????
  if opt.affreg>0  
    VG            = seg8t.tpm(1);
    if ~exist('Ytpmbrain','var')
      Ytpmbrain = spm_read_vols(seg8t.tpm(1)) +  spm_read_vols(seg8t.tpm(2)) +  spm_read_vols(seg8t.tpm(3)); 
    end
    VG.dat(:,:,:) = single(Ytpmbrain); 
    VG.dt         = 16; 
    VG.pinfo      = repmat([1;0],1,size(VG,3));
    VG            = cat_spm_smoothto8bit(VG,6);
  
    VF            = spm_vol(seg8t.image(1));
    VF.dat(:,:,:) = single(Yc{1} + Yc{2} + Yc{3});
    VF.dt         = 16; 
    VF.pinfo      = repmat([1;0],1,size(VF,3));
    VF            = cat_spm_smoothto8bit(VF,6);
    
    if opt.affreg == 2  
      evalc('Affine_com  = cat_vol_set_com(VF);'); % avoid output
      Affine_com(1:3,4) = -Affine_com(1:3,4); %#ok<NODEF> 
    elseif opt.affreg == 3
      Affine_com = eye(4);
    else
      Affine_com = seg8t.Affine; 
    end
  
    % prepare affine parameter 
    aflags = struct('sep',12, ... max(6,max(sqrt(sum(VG(1).mat(1:3,1:3).^2)))), ...
      'regtype','subj','WG',[],'WF',[],'globnorm',1); % job.opts.affreg
    warning off
    Affine  = spm_affreg(VG, VF, aflags, Affine_com); 
    warning on
  elseif opt.affreg<0 || isempty(seg8t.Affine) || all(all(seg8t.Affine==eye(4)))
    %%
    VF            = spm_vol(seg8t.image(1));
    VF.dat(:,:,:) = single(Yc{1} + Yc{2} + Yc{3});
    VF.dt         = 16; 
    VF.pinfo      = repmat([1;0],1,size(VF,3));
    VF            = cat_spm_smoothto8bit(VF,6); %#ok<NASGU> 

    if abs(opt.affreg) == 2 || isempty(seg8t.Affine) || all(all(seg8t.Affine==eye(4))) 
      evalc('Affine_com  = cat_vol_set_com(VF);'); % avoid output
      Affine_com(1:3,4) = -Affine_com(1:3,4); %#ok<NODEF> 
    elseif abs(opt.affreg) == 3
      Affine_com = eye(4);
    else
      Affine_com = seg8t.Affine; 
    end
    warning off
    Affine = spm_maff8(Vo,4,16,seg8t.tpmA,Affine_com,'subj',80); 
    Affine = spm_maff8(Vo,4,4 ,seg8t.tpmA,Affine    ,'subj',40); 
    warning on
  else
    Affine = seg8t.Affine; 
  end
% ##############################
% quantify qc of affine registration (and update)




  % load atlas in individual space by appling the affine transformation
  if ~isempty(opt.Patlas)
    Va = spm_vol(opt.Patlas);
    Ya = zeros(size(Ym),'single'); 
    for zi = 1:size(Ym,3)
      Ya(:,:,zi) = single(spm_slice_vol( Va , ...
        (Va.mat \ Affine * Vo.mat) * spm_matrix([0 0 zi]), ... % apply affine tranformation
        size(Ym,1:2),[0,NaN])); % nearest neighbor interpolation 
    end
    clear Va;
    [~,YD] = cat_vbdist(single(Ya>0),smooth3(Yc{6})<.5); Ya = Ya(YD);
  else
    Ya   = zeros(size(Ym),'single'); 
  end
  
  % load mask in individual space by appling the affine transformation
  if ~isempty(opt.Pmask)
    Vmsk = spm_vol(opt.Pmask);
    Ymsk = zeros(size(Ym),'single'); 
    for zi = 1:size(Ym,3)
      Ymsk(:,:,zi) = single(spm_slice_vol( Vmsk , ...
        (Vmsk.mat \ Affine * Vo.mat) * spm_matrix([0 0 zi]), ... % apply affine tranformation
        size(Ym,1:2),[0,NaN])); % nearest neighbor interpolation 
    end
    clear Vmsk; 
    [~,YD] = cat_vbdist(single(Ymsk>0),smooth3(Yc{6})<.5); Ymsk = Ymsk(YD);
   % Ymsk = Ymsk>1.5; % this mask is limited ... ################## preapare masks ones for faster processing ! #############
  else
    Ymsk = false(size(Ym));
  end
  
  % extend atlas to all voxels
  if ~isempty(opt.Patlas) || ~isempty(opt.Pmask)
    [~,YI] = cat_vbdist(single(Ya>0)); Ya = Ya(YI);   
  end  
  

  % limit boundary box
  Yb  = ( Yc{1} + Yc{2} + Yc{3} ) >.5; 
  [Yo,Ym,Ya,Ymsk,BB] = cat_vol_resize({Yo,Ym,Ya,Ymsk} ,'reduceBrain',tis.res_vx_vol,bd,Yb); 
  for ci = 1:numel(Yc)
    Yc{ci} = cat_vol_resize(Yc{ci} ,'reduceBrain',tis.res_vx_vol,bd,Yb);      
  end

  % limit resolution 
  [Yo,Ym,RES] = cat_vol_resize({Yo,Ym}  ,'reduceV' ,tis.res_vx_vol,opt.reslim,16,'meanm');
  [Ya,Ymsk]   = cat_vol_resize({Ya,Ymsk},'reduceV' ,tis.res_vx_vol,opt.reslim,16,'nearest');
  for ci = 1:numel(Yc)
    Yc{ci} = cat_vol_resize(Yc{ci},'reduceV' ,tis.res_vx_vol,opt.reslim,16,'meanm');
  end
end
%=======================================================================
function [tismri, Ybraindist0]                        = evalSPMseg(Yo,Ym,Yc,Ymsk,vx_vol, fast, opt, seg8t, tis) % TODO: (1) test for WMHs; (2) 
%evalSPMseg. Evaluation of the SPM segmentation.  


  % == evaluate SPM segmentation ==
  % - estimate a brain distance map to further evaluate head/bone/background classes 
  Ybrain         = Yc{1} + Yc{2} + Yc{3}.^2; % underweigt low CSF values 
  [Ybrainr,resR] = cat_vol_resize(Ybrain ,'reduceV' ,vx_vol,3,32,'meanm');  
  Ybraindist0r   = cat_vbdist( single( Ybrainr>.5 ) , true(size(Ybrainr)), resR.vx_volr); 
  Ybraindist0    = cat_vol_resize(Ybraindist0r ,'dereduceV' ,resR);         
  

  % test if the head-class has a lot of very low intensity values
  % the correction is done later if the refinement is selected
  if ~seg8t.isCTseg 
    tismri.warning.clsUpdate = (sum( Yc{5}(:)>.5 & Ym(:)<.3 & Ybraindist0(:)<50 ) / sum( Yc{5}(:)>.5 & Ybraindist0(:)<50 ) ) > .2; 
    if tismri.warning.clsUpdate && 0
      % create a (silent) warning
      cat_io_addwarning( [mfilename ':needSegUpdate'], ...
        sprintf('Bad SPM tissue class 5 - probably overestimated (%0.2f)', ...
        ( sum(Yc{5}(:)>.5 & Ym(:)<.1) ./ sum(Yc{5}(:)>.5 ) ) ) , 1, [1 1], [], 0, opt.verb>1);
    end
  end
% ################### 
% early correction of fatal errors that are easy to correct for ?
% >> no this is refined segmentation 
% ###################

% ###################
% * detect skull-stripping 
% * detect defacing 
% * detect WMHs? (as a second peak in WM class)
% * detect high intesity blood vessels

% ###################

  % separate fat and muscles
  if seg8t.isCTseg % CT images
    Yfat = cat_vol_morph( Yc{5}>0.5 & Ym>-500 & Ym<0   & Ymsk>1 ,'o'); 
    Ymus = cat_vol_morph( Yc{5}>0.5 & Ym>0    & Ym<110 & Ymsk>1 ,'o');
  elseif tis.headFatType==2 % high bone intensity
    Yfat = cat_vol_morph( Yc{5}>0.5 & Ym>mean(tis.seg8n(2))                               & Ymsk>1 ,'o'); 
    Ymus = cat_vol_morph( Yc{5}>0.5 & Ym>mean(tis.seg8n([1,3])) & Ym<mean(tis.seg8n(1:2)) & Ymsk>1 ,'o');
  else % low fat intensity fat suppression
    Yfat = cat_vol_morph( Yc{5}>0.5 & (Ym>mean(tis.seg8n(3))     & Ym<mean(tis.seg8n(1))   |  Ym>mean(tis.seg8n(2))) & Ymsk>1 ,'o'); 
    Ymus = cat_vol_morph( Yc{5}>0.5 &  Ym>mean(tis.seg8n([1,3])) & Ym<mean(tis.seg8n(1:2))                           & Ymsk>1 ,'o');
  end
  


  % == evaluate/test the SPM segmentation == 
  tismri.help.TIV     = 'Total intracranial volume (GM + WM + CSF).';
  tismri.help.vol     = 'Volume of the SPM tissues classes in mm (probability >.5).';
  tismri.help.volr    = 'Relative volume of the SPM tissues classes (probability >.5) normalized by TIV.';
  tismri.help.den     = 'Density of the SPM tissues classes (ie. sum of probability; similar to volume i.e. ~mm).';
  tismri.help.int     = 'Intensity values of the (sub) tissue classes.'; 
  tismri.help.volfat  = 'Volume of fat tissue in the masked upper head (simple threshold to separate the head tissues, in mm).'; 
  tismri.help.volfatr = 'relative volume of fat tissue in the masked upper head (simple threshold to separate the head tissues).'; 
  tismri.help.volmus  = 'Volume of muscle-like tissue in the masked upper head (simple threshold to separate the head tissues, in mm).'; 
  tismri.help.volmusr = 'relative volume of muscle-like tissue in the masked upper head (simple threshold to separate the head tissues).'; 
  tismri.help.denQC   = 'Relation of voxel with high vs. low density within 30 mm distnace. ';
  tismri.help.Tth     = 'Median intensity of the (optimized) tissue class (~peak intensity).'; 
  tismri.help.Tiqr    = 'IQR of the intensity of the (optimized) tissue class (~peak intensity).'; 
  tismri.help.int     = ['Detailed evaluated tissue classes: ' ...
    '(1) GM: median (=tismri.Tth(1)); ' ...
    '(2) WM: median (=tismir.Tth(2)); ' ...
    '(3) CSF: median (=tismri.Th(3)); ' ...
    '(4) bone: kmeans with 2 classes for bone and bone marrow; ' ...
    '(5) head: kmeans with 3 classes head, muscle, fat; ' ...
    '(6) bg: median (=tismri.Th(6)). ' ...
    ];



  tismri.TIV     = sum( (Yc{1}(:) + Yc{2}(:) + Yc{3}(:)) > 0.5) .* prod(vx_vol) / 1000; 
  tismri.volfat  = cat_stat_nansum( Yfat(:) ) .* prod(vx_vol) / 1000; % tissue volume 
  tismri.volmus  = cat_stat_nansum( Ymus(:) ) .* prod(vx_vol) / 1000; % tissue volume 
  tismri.volfatr = tismri.volfat ./ tismri.TIV;
  tismri.volmusr = tismri.volmus ./ tismri.TIV;
  tismri.vol     = nan(1,6); tismri.volr = nan(1,6); tismri.den = nan(1,6); tismri.Tth = nan(1,6);  
  for ci = 1:6
    % estimate tissue volumes and density values
    tismri.vol(ci)  = cat_stat_nansum(Yc{ci}(:)>0.5) .* prod(vx_vol) / 1000; % tissue volume 
    tismri.den(ci)  = cat_stat_nansum(Yc{ci}(:))     .* prod(vx_vol) / 1000; % tissue density
    tismri.volr(ci) = tismri.vol(ci) ./ tismri.TIV;
  
    % test if any class has more low probability values as high
    tismri.cQC(ci) = sum(Yc{ci}(:)>.5) ./ sum(Yc{ci}(:)>eps & Yc{ci}(:)<.5 & Ybraindist0(:)<30); 
    if tismri.cQC(ci)<.5
      cat_io_addwarning( sprintf('%s:badSPMcls%d',mfilename,ci) , ...
        sprintf('Bad SPM tissue class %d - probably underestimated (%0.2f)', ci, tismri.cQC(ci)),1,[1 1],0,0,0);
    end
 
    % refined evaluation with class optimization 
    % otherwise just the median
% ###################    
    switch ci .* (1-fast)
      case 1
        tismri.int.GM   = cat_stat_nanmedian(Yo(Yc{ci}>0.5));
        tismri.Tth(ci)  = tismri.int.GM;
      case 2
        % refinements to avoid PVE, WMHs and PVS  
        Yw                = cat_vol_morph(Yc{ci}>0.9 & cat_vol_morph(Ym>.95,'c'),'e'); 
        if sum(Yw(:)) < 1000, Yw = cat_vol_morph(Yc{ci}>0.5 & cat_vol_morph(Ym>.95,'c'),'e'); end 
        if sum(Yw(:)) < 1000, Yw = cat_vol_morph(Yc{ci}>0.5 & Ym>.9,'e'); end 
        if sum(Yw(:)) < 1000, Yw = Yc{ci}>0.5 & Ym>.9; end 
        tismri.int.WM     = cat_stat_nanmedian(Yo( Yw )); % avoid WMHs and PVE
        tismri.Tth(ci)    = tismri.int.WM;
        % ########################
        % WMHs ? 
        % ########################
        clear Yw; 
      case 3
        % refinements to avoid PVE 
        Yw                = cat_vol_morph(Yc{ci}>0.9 & Ym<.5,'e'); 
        if sum(Yw(:)) < 1000, Yw = Yc{ci}>0.9 | cat_vol_morph(Yc{ci}>0.5 & Ym<.5,'e'); end 
        if sum(Yw(:)) < 1000, Yw = Yc{ci}>0.9 & Ym<.5; end 
        tismri.int.CSF    = cat_stat_nanmedian( Yo( Yw )); % avoid PVE
        tismri.Tth(ci)    = tismri.int.CSF;
        clear Yw;
      case 4
        % ################### need refinement depending on number ###########
        if seg8t.isCTseg 
          tismri.int.bone   = cat_stat_kmeans(Yo(cat_vol_morph(Yc{ci}>0.9,'e')),2); 
          tismri.Tth(ci)    = mean(tismri.int.bone);
        else
          tismri.int.bone   = cat_stat_kmeans(Yo(cat_vol_morph(Yc{ci}>0.9,'e')),3); 
          tismri.Tth(ci)    = tismri.int.bone(3);
        end
       
      case 5
        % ################### need refinement depending on number ###########
        if seg8t.isCTseg 
          tismri.int.head   = cat_stat_kmeans(Yo(cat_vol_morph(Yc{ci}>0.9,'e')),3); 
          tismri.int.muscle = tismri.int.head(2);
          tismri.int.fat    = tismri.int.head(1); 
        else
          tismri.int.head   = cat_stat_kmeans(Yo(cat_vol_morph(Yc{ci}>0.9,'e')),3); 
        end
        tismri.Tth(ci)    = tismri.int.head(3);
      otherwise
        [Yir,Ymr]         = cat_vol_resize( {Yo .* (Yc{ci}>.5),(Yc{ci}>.9)} ,'reduceV' ,1,2,32,'meanm'); 
        tismri.Tth(ci)    = cat_stat_nanmedian(Yir(Ymr>.5));
    end
    % SD in tissue class
    if ~exist('Yir','var')
      [Yir,Ymr]           = cat_vol_resize( {Yo .* (Yc{ci}>.5),(Yc{ci}>.9)} ,'reduceV' ,1,2,32,'meanm'); 
    end
    tismri.Tsd(ci)        = cat_stat_nanstd(Yir(Ymr>.5));
    tismri.Tiqr(ci)       = iqr(Yir(Ymr>.5 & ~isnan(Yir)),'all');
    clear Yir; 
  end
  
  
end
%=======================================================================
function [Yc,Ye,clscor]                               = refineSPM(Yo,Ym,Yc,Ybraindist0,tis,tismri)
% == REFINE SPM DATA ==
% * Create meninges class? 
%   - This is not realy working cause the intensity are changing too much.
%     Hence, it is better to evaluate the whole mid structure. 
%   - Whould thresholding by thickness work? 
%     Expected but it is better to use the thickness as correction factor.
% * What is about blood vessels?
%   - Local filtering ouf outliers might help. 
%     Hence, we use a median filter of bone volume.
% * Split-up skull class by bone and bone marrow?
%   - This was not really working as the conrast changes to much over aging
%     but not robust enough for both structures in various protocols, with 
%     shift artifacts of i.e. fat. 
% * Analyse head by mussles and fat (and skin or other things?)
%   - Not yet.
% * New affine registration for bone only?
%   - should be similar in adults?


% #######################
% Test for errors in SPM segmentation and suggest strategies.  
% * E.g. when the cls4 is also in the real background 
% Strategies: 
% * run SPM with higher seperation (sep)
%
% * create warning in case of to small BB or strong or stange defacing? 

% #######################
% * get FAT peak from SPM to quantify if fat suppression was used
% * get FAT peak from closes head regions 
% >>>  
% * use FAT peak for general scaling 

Yco = Yc; 
vx_vol = tis.res_vx_vol; 
vxmm3  = prod(vx_vol) * 1000; 
Ybraindist0s = cat_vol_smooth3X(Ybraindist0,6);  

  %%
  Yc = Yco; 
  
  % defacing >> background 
  Ydeface = cat_vol_morph(Yo==0,'l',[10,0.1])>0; % need multiple objects for face and ears 
  if sum(Ydeface(:)) > 10000
    for ci = 1:5, Ycn = Yc{ci} .* Ydeface; Yc{6} = Yc{6} + Ycn; Yc{ci} = Yc{ci} - Ycn; end
  end

  % background >> bone | head
  % background should be only outside the head ... not optimal for air
  Yc6    = Yc{6} .* cat_vol_morph(cat_vol_morph(Yc{6}<.5,'lc'),'e'); 
  clscor.BG2BN = sum(Yc6(:) .* ~(Yo(:)<=tismri.Tth(3) & (Ybraindist0s(:)<10))) / vxmm3 / tismri.TIV; 
  clscor.BG2HD = sum(Yc6(:) .*  (Yo(:)<=tismri.Tth(3) & (Ybraindist0s(:)<10))) / vxmm3 / tismri.TIV; 
  Yc{6}  = Yc{6} - Yc6;
  Yc{5}  = Yc{5} + Yc6 .* ~(Yo<=tismri.Tth(3) & (Ybraindist0s<10)); 
  Yc{4}  = Yc{4} + Yc6 .*  (Yo<=tismri.Tth(3) & (Ybraindist0s<10));
  % head >> background 
  % In children the head/bone class is to large and needs correction. 
  if 1 % seg8t.BGtype == ... % low intensity background
    Yc6a = Yc{5} .* smooth3(cat_vol_morph( ((Yc{6} + Yc{5})>.5 & ~Ydeface & Yo<tismri.Tth(3)), 'lo'));
  else
    % noisy background of MP2RAGE/MT sequences need some other (gradient-based) defintion 
  end
  clscor.HD2BN = sum(Yc6a(:)) / vxmm3  / tismri.TIV; 
  Yc{5}  = Yc{5} - Yc6a;
  Yc{6}  = Yc{6} + Yc6a; 
  % final background closeing 
  Yc6    = cat_vol_morph(cat_vol_morph(Yc{6}<.5,'lc')<.5,'e'); 
  for ci = 1:5, Ycn = Yc{ci} .* Yc6; Yc{6} = Yc{6} + Ycn; Yc{ci} = Yc{ci} - Ycn; end


  %% head (+ background) ~ bonemarrow>> bone
  if tis.boneIntType
    Yhead  = min(1,single(Yc{5} + Yc{6} + 0.1 * max(0,Ybraindist0s-15))); % no smoothing here!
    Yhead  = smooth3(cat_vol_morph(cat_vol_morph(smooth3(Yhead)>.7,'lo',4),'d')) .* Yhead .* (Ybraindist0s>0);
    cn = 0; for ci = 5:6, Ycn = Yc{ci} .* Yhead; cn = cn + sum(Ycn(:)); Yc{4} = Yc{4} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
    clscor.HD2BN = cn / vxmm3 / tismri.TIV; 
    Ye{1} = zeros(size(Yc{1}),'single'); 
  else
    dmn = cat_stat_kmeans( Ybraindist0s(Yc{4}>.5 & ~cat_vol_morph(Yc{4}>.5,'e') & Ybraindist0s<(median(Ybraindist0s(Yc{4}>.5))*4) ) , 2 );
    Ybone = cat_vol_morph(cat_vol_morph(Yc{4}>.5 & Ybraindist0s>(dmn(1)*.5) & Ybraindist0s<(dmn(2)*1.2),'lc'),'o',1);
    cn = 0; for ci = [1:3,5:6], Ycn = Yc{ci} .* ~Ybone; cn = cn + sum(Ycn(:)); Yc{4} = Yc{4} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
    clscor.HD2BN = cn / vxmm3 / tismri.TIV; 
   
    %% get the venes
    Yhead  = min(1,single(Yc{5} + Yc{6} + 0.1 * max(0,Ybraindist0s-15))); % no smoothing here!
    Yhead2 = smooth3(cat_vol_morph(cat_vol_morph(smooth3(Yhead)>.7,'lo',4),'d')) .* Yhead .* (Ybraindist0s>0);
    Yhead  = 1 - (Yhead - Yhead2); 
    % Yhead  = cat_vol_morph(Yhead + (Yc{4}>.5),'ldc',3); 
    % Yhead  = cat_vol_smooth3X(Yhead,4)>.5; 
    % Yhead  = Yhead .* cat_vol_morph(cat_vol_morph((Yc{1}+Yc{2}+Yc{3}+Yc{4})>.5,'lc',2),'do',7); 
    %%
    cn = 0; Ye{1} = zeros(size(Yc{1}),'single'); for ci = 5:6, Ycn = Yc{ci} .* Yhead; cn = cn + sum(Ycn(:)); Ye{1} = Ye{1} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
    clscor.HD2BN = cn / vxmm3 / tismri.TIV; 
  end

  %% brain > bone
  Ybrain = single(Yc{1} + Yc{2} + Yc{3}); 
  Ybrain = (cat_vol_morph( Ybrain > .5,'ldo',3)); % & (Ybraindist0>10)) | (Ybrain & (Ybraindist0<10));           % remove skull
  Ybrain = smooth3(cat_vol_morph( Ybrain > 0.5,'dc',3));  % this closing include meninges!
  Ybrain = max( Ybrain , cat_vol_smooth3X(Ybrain,2)>.6); % remove vene
  cn = 0; for ci = 1:3, Ycn = Yc{ci} .* Ybrain; cn = cn + sum(Ycn(:)); Yc{4} = Yc{4} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
  clscor.BR2BN = cn / vxmm3 / tismri.TIV; 
  
  % Yc4 outside the TPM definition >> head | background
  Yc4    = cat_vol_morph(Yc{4}>0,'l') & (Ybraindist0s>10); 
  clscor.BN2HD = sum(Yc4(:)) / vxmm3 / tismri.TIV; 
  Yc{5}  = Yc{5} + Yc{4} .* Yc4; 
  Yc{4}  = Yc{4} - Yc{4} .* Yc4; 

% These correction are not enough in children with inoptimal TPM overlay (eg the NIH templates). 
% In the head worst-case (but good brain case) the bone may include the head fully. 
% Separation maybe by region-growing of (high) head and (low) bone regions. 
% 

  %% head > bone
  %{
  Ybrain = single(Yc{1} + Yc{2} + Yc{3}); 
  Ybrain = (cat_vol_morph( Ybrain > .5,'lo',2) & (Ybraindist0>10)) | (Ybrain & (Ybraindist0<10));           % remove skull
  Ybrain = smooth3(cat_vol_morph( Ybrain > 0.5,'c',3));  % this closing include meninges!
  Ybrain = max( Ybrain , cat_vol_smooth3X(Ybrain,2)>.6); % remove vene
  for ci = 1:3, Ycn = Yc{ci} .* Ybrain; Yc{4} = Yc{4} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
  

  
  %% head >> BG 
  Yc5lth = min(.5,tismri.Tth(5)/tismri.Tth(2));
  Yc5    = Yc{5} .* (Ym<min(0.2,Yc5lth)) .* cat_vol_morph(Yc{5}+Yc{6} > .5,'e'); 
  Yc{6}  = Yc{6} + Yc5; 
  Yc{5}  = Yc{5} - Yc5; 
  %Yc6    = Yc{5} .* (Ym>min(0.2,Yc5lth)) .* cat_vol_morph(Yc{5}+Yc{6},'e'); 
  
  clear Yc5 Yc5lth; 
  %}
  
  %% close bone (by getting from all classes)
  %{
  Yc4    = 1 - (cat_vol_morph(Yc{4}>.25 & ~Ybrain,'lc',5) - Yc{4});
  for ci = [1:3 5:6], Ycn = Yc{ci} .* Yc4; Yc{4} = Yc{4} + (Yc{ci} - Ycn); Yc{ci} = Ycn; end
  %}
% TPM layber model?
% - distance from WM and from BG
%%
  if 0
    %% meninges detection in low int ? 
      Yn  = cat_vol_localstat(single(Ym),(Yc{4}+Yc{3})>.5 & ~cat_vol_morph(Yc{1}>.5,'d'),1,4);
      Yn2 = cat_vol_localstat(single(Ym),(Yc{4}+Yc{3})>.5 & ~cat_vol_morph(Yc{1}>.5,'d'),2,4);
      Ym .* (Yc{4} + Yc{3}),Yn*4 .* Yn2*4 .* (Yc{4}+Yc{3}) .* Ym .* 3; 
      Ym .* (1- smooth3(Yn*4 .* Yn2*4 .* (Yc{4}+Yc{3}).*Ym * 3))
      Ym .* (1- smooth3(Yn*4 .* Yn2*4 .* (Yc{4}+Yc{3}).*Ym * 3)) % ... this is quite nice to remove mengs and BVs (but also affects the GM/WM boudnary !
  end

 
  % ########
  % Problem is to get all bone marrow also the missaligned things that are
  % now part of the head or the brain (aligned as menignes like thing to
  % CSF/GM/WM). 
  %% magic operation that is now working
  test = 0;
  if test
    Ybg = smooth3(Yc{6})>.5 | Ybraindist0s>15 ; 
    Ybe = Ybrain>.5; % & smooth3(Yc{1} + Yc{2} + (Yc{3} .* (Ym<.5)))>.8; % ############## T1! 
    Ybe = cat_vol_morph(Ybe,'lc');
    %Yhd = single( Ybg ); Yhd( Ybe ) = nan; [~,Yhd] = cat_vol_downcut( Yhd , max(eps,min(10,  3 - Ym*3)), .5);
    %Yhd = single( Ybg ); Yhd( Ybe ) = nan; [~,Yhd] = cat_vol_downcut( Yhd , max(eps,min(10,  Ym*3)), .5);
    Ybd = single( Ybe ); Ybd( Ybg ) = nan; [~,Ybd] = cat_vol_downcut( Ybd , min(10,  Ym*3) , .5);
    %Yboneh = ( (min(Yhd,Ybd*1e10) ./ (Yhd+Ybd))>.5);
    Yboneh = Ybd./max(1,Ybraindist0s - median(Ybraindist0s(Yc{4}(:)>.5))) / 100; 
    Yboneh(Ybg | isnan(Yboneh) | Yboneh>2) = 0;  
    Ych = zeros(Yboneh);
    for ci = [1:3 5:6], Ych = Ych + Yc{5} .* Yboneh; end
    Yc5   = Yc{5} .* Yboneh;
    Yc{4} = Yc{4} + Yc5;
    Yc{5} = Yc{5} - Yc5; 
    clear Yc5; 
  end

  
end
%=======================================================================
function [Ybonepp, Ybonethick, Ybonemarrow, Yheadthick, val] = extractbone(Vo,Yo,Ym,Yc,Ye,Ya,Ymsk,seg8t,tis,out,opt,vx_vol,RES,BB)
%% * Report: 
%   - better an upper slice?
%   - optimize print font size
%   - add basic parameters (reduce,mask,Atlas,tpm)
%   - add tissue volumes
%   - add tissue intensities
%   - add histograms for tissues
%   - modify colorbar position and labeling > CAT
%   - use fat/musle colors (yellow, pink)  
%   - use green/cyan for bone?
%   - affine registration surf problem
%   - use final segmenation for overlay but mark outliers 
%   - Opimize report line >> table like similar to QC with vols, thick & intensities 

  %%
  if tis.weighting == -1
  %% CT images
    Ybrain0      = cat_vol_morph(cat_vol_morph((Yc{1} + Yc{2} + Yc{3})>.5,'lc',1,vx_vol),'lo',3,vx_vol); % remove 
    Ybraindist1  = cat_vbdist( single(Ybrain0>0.5) , (Yc{4} + Yc{5})>0 , vx_vol);
    Yhead        = Yc{1} + Yc{2} + Yc{3} + Yc{4}; Yhead(cat_vol_morph(Yhead>.5,'ldc',2)) = 1; 
    Ybone        = Yc{4};
    Ybrain       = (Yhead - Ybone) .* cat_vol_morph(Yhead>.5,'e') .* (Ybraindist1<4); % .* Ybrain;
    Ybraindist   = cat_vbdist( single(Ybrain>0.5) , Ybone>.5, vx_vol) .* (Ybone>0.5);
    Yheaddist    = cat_vbdist( single(Yhead<0.5)  , Ybone>.5, vx_vol) .* (Ybone>0.5);
    Ybonethick   = Ybraindist  + Yheaddist;  % correct for voxel-size
    Ybonepp      = min(1,Yheaddist  ./ max(eps,Ybonethick));  Ybonepp(Ybrain>.5) = 1; % percentage map to
    Ybonemarrow  = Ym;

  else
  % MRI iamges
    if 0
      Yheadbone  = cat_vol_morph(cat_vol_morph( smooth3(Yc{4} + Yc{5} + Yc{6}) > 0.5,'c',3),'e');     % create 
      Yheadbone  = cat_vol_smooth3X(Yheadbone,4)>.5;
      Ybrain     = single(Yc{1} + Yc{2} + Yc{3} - Yheadbone);                                 % SPM brain 
      Ybrainc    = (Yc{1} + Yc{2} + Yc{3}) .* (Ybrain<.5 & Yheadbone); 
      Yc{4}      = Yc{4} + Ybrainc; 
      Ybrain     = max( Ybrain , cat_vol_smooth3X(Ybrain,4)>.5); % remove vene
      %%
      Yhead      = min(1,single(Yc{5} + Yc{6}));
      Yhead(~cat_vol_morph(cat_vol_morph(Yhead>.5,'lo',1),'d') & Yhead>.5 & Ybrain<.5) = 0.5; 
      Ybg        = single(Yc{6});
      Ybone      = single(Yc{4});
     
      %% 
      if 1
        Ybrain = single(Ybrain); Ybone = single(Ybone); Yhead = single(Yhead);
        spm_smooth(Ybrain,Ybrain,2./vx_vol);
        spm_smooth(Yhead,Yhead,2./vx_vol);
        spm_smooth(Ybg,Ybg,2./vx_vol);
      end
      Ybone     = ~(Ybrain>.5 | Yhead>.5); 
    else
      Ybrain     = single(Yc{1} + Yc{2} + Yc{3});
      Yhead      = single(Yc{5} + Yc{6});
      Ybone      = single(Yc{4});
    end
    if ~isempty(Ye)
      Ybrain = Ybrain + Ye{1};
    end
    
    % remove PVE
    %  Ybonee     = Ybone>.5 & ~cat_vol_morph(Ybone>.5,'e');
    %  Ybonem     = cat_vol_localstat(single(Yo + ((Ybone<.5)*100)),Ybone>.5,2,2); 
    
    Ybrain = cat_vol_morph( cat_vol_morph( Ybrain>.5 , 'lc' , 2) , 'lo' , 2);
    Yhead  = ~cat_vol_morph( Yhead<.5 , 'lo' , 1); 

    Ybrain = cat_vol_smooth3X( Ybrain , 1)>.5; 
    Yhead  = cat_vol_smooth3X( Yhead  , 1)>.5; 
    Ybone1 = 1 - Ybrain - Yhead; 
    

    %% bone layers
    Ybraindist   = cat_vbdist( single(Ybrain>0.5) , Ybone1>0, vx_vol);
    Yheaddist    = cat_vbdist( single(Yhead>0.5)  , Ybone1>0, vx_vol);
    Ybonethick   = Ybraindist  + Yheaddist;  % correct for voxel-size
    Ybonepp      = min(1,Yheaddist  ./ max(eps,Ybonethick)); Ybonepp(Ybrain>.5) = 1; % percentage map to
    if 0
      % head values
      Ybonehead = ~(Ybrain>.5 | Ybg>.5); 
      %Ybonehead    = Ybonehead .* (Ybraindist<30 & Ybonehead); 
      % Ybg          = max(Ybg,Ybraindist>=30 & Ybonehead); 
      Ybraindist2  = cat_vbdist( single(Ybrain>0.5) , Ybonehead, vx_vol) .* (Ybonehead>0.5);
      Yheaddist2   = cat_vbdist( single(Ybg>0.5)    , Ybonehead, vx_vol) .* (Ybonehead>0.5);
      Ybonethick2  = Ybraindist2 + Yheaddist2;  % correct for voxel-size
      Ybonepp2     = min(1,Yheaddist2 ./ max(eps,Ybonethick2)); Ybonepp2(Ybrain>.5) = 1;    % percentage map to
      Ybonepp      = max(Ybonepp,max(0,Ybonepp2*3-2));
    end
    clear braindist Ybonedist
    
    % bonemarrow
    Ybonemarrow = single(Ym/tis.seg8n(3)) .* (Ybone>.5); % bias intensity normalized corrected
  end 


  % head 
  Yskull       = single(Ym/tis.seg8n(3)) .* (Yc{5}>.5);
  Ybgdist      = cat_vbdist( single(Yc{6}) , Ybrain<0.5, vx_vol);
  Ybndist      = cat_vbdist( cat_vol_morph(single(Ybrain + Ybone),'lc') , Yc{6}<.5, vx_vol);
  Yheadthick   = Ybndist + Ybgdist - max(0,Yheaddist - .5);
  Yheadthick   = cat_vol_localstat(Yheadthick,Yheadthick<1000,1,1);
  [~,YD]       = cat_vbdist(single(Yheadthick>0)); Yheadthick = Yheadthick(YD); % simple extension to avoid NaNs (BB limit is ok)
  clear Ybgdist Ybndist


 

%% ###################
% edge-measure ! 
% * this one is not realy working (in low fat cases?)
% * bone / bone-marrow segmenation: 
%   - bone marrow as outstanding maximum structure in the middle/center area
%   


  %% measures as column elements
  rii = 1;
  val.help = 'ROI=0 is defined masked global values excluding the lower parts of the skull, whereas all other ROIs are without masking';
  for ri = 0:max(Ya(Ya(:)<intmax('uint16')))
    if ri == 0 || isnan(ri)
      ri = 0; %#ok<FXSET> % case of failed atlas mapping 
      val.boneatlas_id(1,rii)       = inf;
      val.nonnanvol(1,rii)          = sum(Ya(:)>intmax('uint16')) ./ numel(Ya(:));
      if opt.mask, val.boneatlas_name{1,rii} = 'full-masked'; 
      else,        val.boneatlas_name{1,rii} = 'full-unmasked'; 
      end
      % bone (marrow) intensity >> rename later to skull (skull = bone + bone-marrow)  
      val.bonemarrow_mean(1,rii)    = cat_stat_nanmean(   Ybonemarrow( Ymsk(:)>1 & Ybonemarrow(:)~=0 ) ); 
      val.bonemarrow_std(1,rii)     = cat_stat_nanstd(    Ybonemarrow( Ymsk(:)>1 & Ybonemarrow(:)~=0 ) ); 
      val.bonemarrow_med(1,rii)     = cat_stat_nanmedian( Ybonemarrow( Ymsk(:)>1 & Ybonemarrow(:)~=0 ) ); 
      val.bonemarrow_iqr(1,rii)     = iqr(                Ybonemarrow( Ymsk(:)>1 & Ybonemarrow(:)~=0 ) ); 
      % bone thickness
      val.bonethickness_mean(1,rii) = cat_stat_nanmean(   Ybonethick( Ymsk(:)>1  & Ybonethick(:)~=0 ) ); 
      val.bonethickness_std(1,rii)  = cat_stat_nanstd(    Ybonethick( Ymsk(:)>1  & Ybonethick(:)~=0 ) ); 
      val.bonethickness_med(1,rii)  = cat_stat_nanmedian( Ybonethick( Ymsk(:)>1  & Ybonethick(:)~=0 ) ); 
      val.bonethickness_iqr(1,rii)  = iqr(                Ybonethick( Ymsk(:)>1  & Ybonethick(:)~=0 ) ); 
      % head thickness
      val.head_mean(1,rii)          = cat_stat_nanmean(   Yskull( Ymsk(:)>1  & Yskull(:)~=0 ) ); 
      val.head_std(1,rii)           = cat_stat_nanstd(    Yskull( Ymsk(:)>1  & Yskull(:)~=0 ) ); 
      val.head_med(1,rii)           = cat_stat_nanmedian( Yskull( Ymsk(:)>1  & Yskull(:)~=0 ) ); 
      val.head_iqr(1,rii)           = iqr(                Yskull( Ymsk(:)>1  & Yskull(:)~=0 ) ); 
      % head thickness
      val.headthickness_mean(1,rii) = cat_stat_nanmean(   Yheadthick( Ymsk(:)>1  & Yheadthick(:)~=0 ) ); 
      val.headthickness_std(1,rii)  = cat_stat_nanstd(    Yheadthick( Ymsk(:)>1  & Yheadthick(:)~=0 ) ); 
      val.headthickness_med(1,rii)  = cat_stat_nanmedian( Yheadthick( Ymsk(:)>1  & Yheadthick(:)~=0 ) ); 
      val.headthickness_iqr(1,rii)  = iqr(                Yheadthick( Ymsk(:)>1  & Yheadthick(:)~=0 ) ); 
      rii = rii + 1;
    else
      if sum(Ya(:)==ri)~=0
        val.boneatlas_id(1,rii)       = ri;  
        val.boneatlas_name{1,rii}     = sprintf('ROI%d',ri); 
        val.nonnanvol(1,rii)          = sum(Ya(:)==ri) ./ numel(Ya(:));
        % bone marrow intensity 
        val.bonemarrow_mean(1,rii)    = cat_stat_nanmean(   Ybonemarrow( Ymsk(:)>1 & Ybonemarrow(:)~=0 & Ya(:)==ri) ); 
        val.bonemarrow_std(1,rii)     = cat_stat_nanstd(    Ybonemarrow( Ymsk(:)>1 & Ybonemarrow(:)~=0 & Ya(:)==ri) );
        val.bonemarrow_med(1,rii)     = cat_stat_nanmedian( Ybonemarrow( Ymsk(:)>1 & Ybonemarrow(:)~=0 & Ya(:)==ri) );
        val.bonemarrow_iqr(1,rii)     = iqr(                Ybonemarrow( Ymsk(:)>1 & Ybonemarrow(:)~=0 & Ya(:)==ri) );
        % thickness
        val.bonethickness_mean(1,rii) = cat_stat_nanmean(   Ybonethick( Ymsk(:)>1  & Ybonethick(:)~=0  & Ya(:)==ri) );
        val.bonethickness_std(1,rii)  = cat_stat_nanstd(    Ybonethick( Ymsk(:)>1  & Ybonethick(:)~=0  & Ya(:)==ri) );
        val.bonethickness_med(1,rii)  = cat_stat_nanmedian( Ybonethick( Ymsk(:)>1  & Ybonethick(:)~=0  & Ya(:)==ri) );
        val.bonethickness_iqr(1,rii)  = iqr(                Ybonethick( Ymsk(:)>1  & Ybonethick(:)~=0  & Ya(:)==ri) );
        % head intensity
        val.head_mean(1,rii)          = cat_stat_nanmean(   Yskull( Ymsk(:)>1  & Yskull(:)~=0  & Ya(:)==ri) ); 
        val.head_std(1,rii)           = cat_stat_nanstd(    Yskull( Ymsk(:)>1  & Yskull(:)~=0  & Ya(:)==ri) ); 
        val.head_med(1,rii)           = cat_stat_nanmedian( Yskull( Ymsk(:)>1  & Yskull(:)~=0  & Ya(:)==ri) );
        val.head_iqr(1,rii)           = iqr(                Yskull( Ymsk(:)>1  & Yskull(:)~=0  & Ya(:)==ri) );
        % head thickness
        val.headthickness_mean(1,rii) = cat_stat_nanmean(   Yheadthick( Ymsk(:)>1  & Yheadthick(:)~=0  & Ya(:)==ri) ); 
        val.headthickness_std(1,rii)  = cat_stat_nanstd(    Yheadthick( Ymsk(:)>1  & Yheadthick(:)~=0  & Ya(:)==ri) ); 
        val.headthickness_med(1,rii)  = cat_stat_nanmedian( Yheadthick( Ymsk(:)>1  & Yheadthick(:)~=0  & Ya(:)==ri) );
        val.headthickness_iqr(1,rii)  = iqr(                Yheadthick( Ymsk(:)>1  & Yheadthick(:)~=0  & Ya(:)==ri) );
        % addroi
        rii = rii + 1;
      end
    end
  end





  %% restore resolution & boundary box
  [Ybonepp, Ybonethick, Ybonemarrow,Yheadthick] = cat_vol_resize({Ybonepp, Ybonethick, Ybonemarrow,Yheadthick} ,'dereduceV' ,RES); % ############### INTERPOLATION ???
  [Ybonepp, Ybonethick, Ybonemarrow,Yheadthick] = cat_vol_resize({Ybonepp, Ybonethick, Ybonemarrow,Yheadthick} ,'dereduceBrain',BB); 

  if tis.boneIntType == 0 && tis.weighting > 0
    Ybonemarrow = Ybonemarrow * 3;
  end
  if tis.weighting == -1
    Ybonemarrow(Ybonemarrow==0) = -1024;
  end

  


  %% write output maps
  %  - what maps do we really need?
  %  - intensity normalized maps used for normalization and analysis? 
  %  - 
  if opt.writevol
    
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
    job.output.position    = struct('native',0,'warped',0,'dartel',0);
    job.output.bone        = struct('native',0,'warped',0,'dartel',0);
    job.output.bonethick   = struct('native',0,'warped',0,'dartel',0);
    job.output.headthick   = struct('native',0,'warped',0,'dartel',0);

    % midline map also for masking masking
    %cat_io_writenii(Vo,Ybonemid,out.P.mridir,sprintf('bonemid%d_',opt.method), ...
    %  'bone percentage position midline map','uint8',[0,1/255], ... 
    %  min([1 0 2],[job.output.position.native job.output.position.warped job.output.position.dartel]),trans); 
    % masked map for averaging
    cat_io_writenii(Vo,Ybonepp,out.P.mridir,sprintf('bonepp%d_',opt.method), ...
      'bone percentage position map','uint16',[0,0.001], ... 
      min([1 0 2],[job.output.bone.native job.output.bone.warped job.output.bone.dartel]),trans);
    cat_io_writenii(Vo,Ybonemarrow,out.P.mridir,sprintf('bonemarrow%d_',opt.method), ...
      'bone percentage position map','uint16',[0,0.001], ... 
      min([1 0 2],[job.output.bonemarrow.native job.output.bonemarrow.warped job.output.bonemarrow.dartel]),trans);
    cat_io_writenii(Vo,Ybonethick,out.P.mridir,sprintf('bonethickness%d_',opt.method), ...
      'bone thickness map','uint16',[0,0.001], ... 
      min([1 0 2],[job.output.bonethick.native job.output.bonethick.warped job.output.bonethick.dartel]),trans);
    cat_io_writenii(Vo,Yheadthick,out.P.mridir,sprintf('headthickness%d_',opt.method), ...
      'head thickness map','uint16',[0,0.001], ... 
      min([1 0 2],[job.output.headthick.native job.output.headthick.warped job.output.headthick.dartel]),trans);
  end
end
%=======================================================================
function [Si, St, Stm, Sth, sROI]                                = create_bone_surface(Vo,Ym,Yc,Ybonepp,Ybonethick,Ybonemarrow,Yheadthick,Ya,Ymsk,out,opt)
%create_bone_surface. Surface-bases processing pipeline.
  vx_vol = sqrt(sum(Vo.mat(1:3,1:3).^2));
  
  % surface coordinate transformation matrix
  matlab_mm  = Vo.mat * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];   % CAT internal space
 

  %% create surface 
  Ypp  = Ybonepp + 0; spm_smooth(Ypp,Ypp,4);  %Ypp  = Yc{5};  spm_smooth(Ypp,Ypp,4);
  [Yboneppr,res] = cat_vol_resize(smooth3(Ybonepp),'reduceV',vx_vol,opt.reduce,6,'meanm'); %#ok<ASGLU>        ./min(vx_vol)
  txt = evalc(sprintf('[Yppc,CBS.faces,CBS.vertices] = cat_vol_genus0(Yboneppr,.5,0);')); %#ok<NASGU>  % try to use lower values to avoid running into the vene
  CBS.vertices = CBS.vertices .* repmat(res.vx_red,size(CBS.vertices,1),1) - repmat((res.vx_red-1)/2,size(CBS.vertices,1),1); %#ok<NODEF> 
  % CBS.vertices = CBS.vertices .* res.vx_red(1) - (res.vx_red-1)/2; %#ok<NODEF> 
  CBS = cat_surf_fun('smat',CBS,matlab_mm); % transform to mm 
  CBS.EC = size(CBS.vertices,1) + size(CBS.faces,1) - size(spm_mesh_edges(CBS),1);  
  saveSurf(CBS,out.P.central);
  % optimize surface for midbone position 
  % simple blurring 
  Vpp  = cat_io_writenii(Vo, Ypp , '', sprintf('%s.pp','skull') ,...
      'percentage bone position map', 'uint8', [0,1/255],[1 0 0],struct());
  cmd = sprintf('CAT_BlurSurfHK "%s" "%s" %d', out.P.central ,out.P.central, 16 );  % simple smoothing to remove stair artifacts - 8-16 iterations in red 2
  cat_system(cmd,0);
  if 0 % not working?
    cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...   
                   'avg  %0.3f  %0.3f .2  .1  5  0 "0.5"  "0.5"  n 0  0  0 %d %g  0.0 0'], ...          
                    Vpp.fname,out.P.central ,out.P.central ,-0.00005, 0.00005, 100, 0.01); %#ok<UNRCH> % -0.05, 0.05
    cat_system(cmd,0);
  end
  CBS = loadSurf(out.P.central);
  if exist(Vpp.fname,'file'), delete(Vpp.fname); end


%% extract values from the bone 
  % create a smoothed thickness map for the mapping
  Ybonethick2 = cat_vol_approx(Ybonethick .* (Ybonethick>1 & Ybonethick<100),'nh',1,4);
  Ybonethick(Ybonethick>100) = Ybonethick2(Ybonethick>100);
  spm_smooth(Ybonethick,Ybonethick2,2./vx_vol);
  % write smoothed thickness map
  Si = CBS; Si.facevertexcdata = cat_surf_fun('isocolors',max(3,Ybonethick), CBS, matlab_mm); 
  M  = spm_mesh_smooth(Si);
  % add min
  Si.facevertexcdata = single(spm_mesh_smooth(M,double(Si.facevertexcdata),10));
  Si.facevertexcdata(Si.facevertexcdata<=1) = nan; 
  Si = cat_surf_fun('approxnans',Si);

  cat_io_FreeSurfer('write_surf_data',out.P.thick,Si.facevertexcdata);

  %% map values
  % write bone marrow map 
  Vpp  = cat_io_writenii(Vo, Ybonemarrow , '', 'skull.marrow' ,'bone marrow', 'single', [0,1],[1 0 0],struct());
  estminmax = 1; 
  if estminmax
    %% estimate the minimum value > real bone
    Ybonemarrow3 = Ybonemarrow; Ybonemarrow3(Ybonepp==0 | Ybonepp==1) = cat_stat_nanmean(Ybonemarrow3(Ybonepp>0 & Ybonepp<1)); 
    Vppmin  = cat_io_writenii(Vo, Ybonemarrow3 , '', 'skull.bone' ,'bone', 'single', [0,1],[1 0 0],struct());
    mappingstr = sprintf('-linear -min -steps "10" -start "-.5" -end ".5" -thickness "%s" ', out.P.thick); % weighted_avg
    cmd = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s" ',mappingstr, out.P.central,  Vppmin.fname , out.P.marrowmin );
    cat_system(cmd,0); delete(Vppmin.fname);
    %% estimate the maximum value > bone marrow
    mappingstr = sprintf('-linear -max -steps "10" -start "-0.25" -end "0.25" -thickness "%s" ', out.P.thick); % weighted_avg
    cmd = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s" ',mappingstr, out.P.central,  Vpp.fname , out.P.marrowmax );
    cat_system(cmd,0);
    marrowmin = cat_io_FreeSurfer('read_surf_data',out.P.marrowmin);
    marrowmax = cat_io_FreeSurfer('read_surf_data',out.P.marrowmax);
  end
  %%
  mappingstr = sprintf('-linear -weighted_avg -steps "10" -start "-.25" -end ".25" -thickness "%s" ', out.P.thick); % weighted_avg
  cmd = sprintf('CAT_3dVol2Surf %s "%s" "%s" "%s" ',mappingstr, out.P.central,  Vpp.fname , out.P.marrow );
  cat_system(cmd,0);
  %%
  Si.facevertexcdata = cat_io_FreeSurfer('read_surf_data',out.P.marrow);
  % get atlas information
  Satlas = cat_surf_fun('isocolors',Ya, CBS, matlab_mm,'nearest');
  Si.facevertexcdata = Si.facevertexcdata .* (Satlas>0);
  cat_io_FreeSurfer('write_surf_data',out.P.thick,Si.facevertexcdata);

  if 0
    %%
    Si.facevertexcdata = cat_io_FreeSurfer('read_surf_data',out.P.marrowmax); %#ok<UNRCH> 
    %Si.facevertexcdata = single(spm_mesh_smooth(M,double(Si.facevertexcdata),2));
    %facevertexcdata = single(spm_mesh_smooth(M,double(cat_io_FreeSurfer('read_surf_data',out.P.marrowmin)),2));
    Si.facevertexcdata = Si.facevertexcdata .* double(cat_io_FreeSurfer('read_surf_data',out.P.marrowmin));
    Sht = cat_surf_render2(Si); title(sprintf('bone thickness %s/%s',spm_str_manip(pp,'t'),ff));
    cat_surf_render2('Clim',Sht,[0 8]); cat_surf_render2('Colorbar',Sht);
  end
  %%
  if 0
    %% fat
    Yppi = max(0,min(1, smooth3(cat_vol_morph(Ybone<.5,'e') + min(0,-Ydiv*4)) + smooth3(Ybonepp))); %#ok<UNRCH> 
    Vpp  = cat_io_writenii(Vo, Yppi , '', sprintf('%s.pp','skull') ,...
        'percentage position map - V2 - pial=1/3, central=1/2, white=2/3, 0 and 1 to stablize white/pial values', 'uint8', [0,1/255],...
        min([1 1 2],[1 0 0]),struct());
    cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...   
                   'avg  %0.3f  %0.3f .2  .1  5  0 "0.1"  "0.1"  n 0  0  0 %d %g  0.0 0'], ...          
                    Vpp.fname,Pouter,Pouter,-.1, .1, 50, 0.01); 
    cat_system(cmd,opt.verb-3);
    CBSo = loadSurf(Pouter);
  end








  %% get peak bone threshold
  % -use of nan for unwanted rois ?
  Si.vertices = CBS.vertices; Si.faces = single(CBS.faces); St = Si; Stm = Si; Sth = Si; Smsk = Si; 
  S.atlas     = cat_surf_fun('isocolors',Ya         , CBS, matlab_mm, 'nearest');
  S.thick     = cat_surf_fun('isocolors',Ybonethick , CBS, matlab_mm);
  S.hdthick   = cat_surf_fun('isocolors',Yheadthick , CBS, matlab_mm);
  St.facevertexcdata = S.thick;
  St.facevertexcdata(St.facevertexcdata<=2) = nan; 
  Sth.facevertexcdata = S.hdthick; 
  Sth.facevertexcdata(Sth.facevertexcdata<=2) = nan; 
  St  = cat_surf_fun('approxnans',St);
  Si  = cat_surf_fun('approxnans',Si);
  Stm.facevertexcdata = S.thick .* cat_surf_fun('isocolors',max(.1,1 - (cat_vol_grad(Ya*1000)>0.1) * .9 ),CBS, matlab_mm); 
  Smsk.facevertexcdata( cat_surf_fun('isocolors', single( Ymsk ) ,CBS, matlab_mm) > 1.5) = nan; 
  Sth.facevertexcdata( isnan( Smsk.facevertexcdata) ) = nan; 
  if opt.mask
    mskROIs = 0; % define unwanted ROIs
    for mski = mskROIs, Si.facevertexcdata(S.atlas==mski) = nan; end
  end
  cat_io_FreeSurfer('write_surf_data',out.P.thick,St.facevertexcdata);
  cat_io_FreeSurfer('write_surf_data',out.P.headthick,St.facevertexcdata);
  if 0
    %% just to test surfaces and values here
    Sx.vertices = CBS.vertices; Sx.faces = single(CBS.faces); 
    switch 6 %#ok<UNRCH> 
      case 0, Sx.facevertexcdata = cat_surf_fun('isocolors',Ybonemarrow, CBS, matlab_mm); 
      case 1, Sx.facevertexcdata = cat_surf_fun('isocolors',Ybonethick,CBS, matlab_mm); 
      case 2, Sx.facevertexcdata = cat_surf_fun('isocolors',Ym .* max(0, 1 + (cat_vol_grad(Ya*1000)>0.1) * 10 ),CBS, matlab_mm); 
      case 3, Sx.facevertexcdata = cat_surf_fun('isocolors',Ya,CBS, matlab_mm,'nearest'); 
      case 4, Sx.facevertexcdata = cat_surf_fun('isocolors',Ybgdist - Yheaddist,CBS, matlab_mm,'nearest'); 
      case 5, Sx.facevertexcdata = cat_surf_fun('isocolors',Yheadthick,CBS, matlab_mm); 
      case 6, Sx.facevertexcdata = Sth.facevertexcdata; 
    end  

    Sht = cat_surf_render2(Sx); title(sprintf('bone thickness %s/%s',strrep(spm_str_manip(out.P.pp,'t'),'_','\_'),out.P.orgff));
    cat_surf_render2('Clim',Sht,[0 4]); cat_surf_render2('Colorbar',Sht);
  end

  
  %% measures as column elements
  rii = 1;
  sROI.help = 'ROI=0 is defined masked global values excluding the lower parts of the skull, whereas all other ROIs are without masking';
  for ri = 0:max(Ya(Ya(:)<intmax('uint16')))
    if ri == 0 || isnan(ri)
      ri = 0; %#ok<FXSET> % case of failed atlas mapping 
      sROI.boneatlas_id(1,rii)       = inf;
      sROI.nonnanvol(1,rii)          = sum(S.atlas>intmax('uint16')) ./ numel(S.atlas);
      if opt.mask, sROI.boneatlas_name{1,rii} = 'full-masked'; 
      else,        sROI.boneatlas_name{1,rii} = 'full-unmasked'; 
      end
      % bone marrow intensity 
      sROI.bonemarrow_mean(1,rii)    = cat_stat_nanmean(Si.facevertexcdata .* Smsk.facevertexcdata); 
      sROI.bonemarrow_std(1,rii)     = cat_stat_nanstd(Si.facevertexcdata .* Smsk.facevertexcdata); 
      sROI.bonemarrow_med(1,rii)     = cat_stat_nanmedian(Si.facevertexcdata .* Smsk.facevertexcdata); 
      sROI.bonemarrow_iqr(1,rii)     = iqr(Si.facevertexcdata .* Smsk.facevertexcdata); 
      % bone thickness
      sROI.bonethickness_mean(1,rii) = cat_stat_nanmean(S.thick); 
      sROI.bonethickness_std(1,rii)  = cat_stat_nanstd(S.thick); 
      sROI.bonethickness_med(1,rii)  = cat_stat_nanmedian(S.thick); 
      sROI.bonethickness_iqr(1,rii)  = iqr(S.thick); 
      % head thickness
      sROI.headthickness_mean(1,rii) = cat_stat_nanmean(S.hdthick); 
      sROI.headthickness_std(1,rii)  = cat_stat_nanstd(S.hdthick); 
      sROI.headthickness_med(1,rii)  = cat_stat_nanmedian(S.hdthick); 
      sROI.headthickness_iqr(1,rii)  = iqr(S.thick); 
      if estminmax
        sROI.bonemarrowmin_mean(1,rii)    = cat_stat_nanmean(marrowmin); 
        sROI.bonemarrowmin_std(1,rii)     = cat_stat_nanstd(marrowmin); 
        sROI.bonemarrowmin_med(1,rii)     = cat_stat_nanmedian(marrowmin); 
        sROI.bonemarrowmin_iqr(1,rii)     = iqr(marrowmin); 
        sROI.bonemarrowmax_mean(1,rii)    = cat_stat_nanmean(marrowmax); 
        sROI.bonemarrowmax_std(1,rii)     = cat_stat_nanstd(marrowmax); 
        sROI.bonemarrowmax_med(1,rii)     = cat_stat_nanmedian(marrowmax); 
        sROI.bonemarrowmax_iqr(1,rii)     = iqr(marrowmax); 
      end
      rii = rii + 1;
    else
      if sum(S.atlas==ri)>0
        sROI.boneatlas_id(1,rii)       = ri;  
        sROI.boneatlas_name{1,rii}     = sprintf('ROI%d',ri); 
        sROI.nonnanvol(1,rii)          = sum(S.atlas==ri) ./ numel(S.atlas);
        % bone marrow intensity 
        sROI.bonemarrow_mean(1,rii)    = cat_stat_nanmean(Si.facevertexcdata(S.atlas==ri)); 
        sROI.bonemarrow_std(1,rii)     = cat_stat_nanstd(Si.facevertexcdata(S.atlas==ri)); 
        sROI.bonemarrow_med(1,rii)     = cat_stat_nanmedian(Si.facevertexcdata(S.atlas==ri)); 
        sROI.bonemarrow_iqr(1,rii)     = iqr(Si.facevertexcdata(S.atlas==ri)); 
        % thickness
        sROI.bonethickness_mean(1,rii) = cat_stat_nanmean(S.thick(S.atlas==ri)); 
        sROI.bonethickness_std(1,rii)  = cat_stat_nanstd(S.thick(S.atlas==ri)); 
        sROI.bonethickness_med(1,rii)  = cat_stat_nanmedian(S.thick(S.atlas==ri)); 
        sROI.bonethickness_iqr(1,rii)  = iqr(S.thick(S.atlas==ri)); 
        % head thickness
        sROI.headthickness_mean(1,rii) = cat_stat_nanmean(S.hdthick(S.atlas==ri)); 
        sROI.headthickness_std(1,rii)  = cat_stat_nanstd(S.hdthick(S.atlas==ri)); 
        sROI.headthickness_med(1,rii)  = cat_stat_nanmedian(S.hdthick(S.atlas==ri)); 
        sROI.headthickness_iqr(1,rii)  = iqr(S.thick(S.atlas==ri)); 
        if estminmax
          sROI.bonemarrowmin_mean(1,rii)    = cat_stat_nanmean(marrowmin(S.atlas==ri)); 
          sROI.bonemarrowmin_std(1,rii)     = cat_stat_nanstd(marrowmin(S.atlas==ri)); 
          sROI.bonemarrowmin_med(1,rii)     = cat_stat_nanmedian(marrowmin(S.atlas==ri)); 
          sROI.bonemarrowmin_iqr(1,rii)     = iqr(marrowmin(S.atlas==ri)); 
          sROI.bonemarrowmax_mean(1,rii)    = cat_stat_nanmean(marrowmax(S.atlas==ri)); 
          sROI.bonemarrowmax_std(1,rii)     = cat_stat_nanstd(marrowmax(S.atlas==ri)); 
          sROI.bonemarrowmax_med(1,rii)     = cat_stat_nanmedian(marrowmax(S.atlas==ri)); 
          sROI.bonemarrowmax_iqr(1,rii)     = iqr(marrowmax(S.atlas==ri)); 
        end
        rii = rii + 1;
      end
    end
  end
end
%=======================================================================
function out = simpleBone(seg8t,Ym,Yp)
%simpleBone. First simple version that worked quite well.

% == (1) SPM gaussians ==
  % get high prob value of CSF tissue
  out.spm_con = sum(seg8t.mn( seg8t.lkp == 3 )' .* ...
    ( seg8t.mg( seg8t.lkp == 3 ) == max(seg8t.mg( seg8t.lkp == 3 ) ) ));
  
  if max(seg8t.lkp) > 6
    out.spm_bone_med1 = seg8t.mn( seg8t.lkp == 4 ) / out.spm_con;
    out.spm_bone_med2 = seg8t.mn( seg8t.lkp == 5 ) / out.spm_con;
  else
    % default case for 6 TPM classes
    bone = sort(seg8t.mn( seg8t.lkp == 4 ) / out.spm_con); % order intensities within class 4 (normalized by CSF)
    out.spm_bone_med1 = bone(min(numel(bone),3)); % highest value
    out.spm_bone_med2 = bone(1); 
    out.spm_bone_med3 = bone(min(numel(bone),2)); 
  end
 

% == (2) median bone intensity ==
  % load SPM tissue segments
  for ci = 1:5 
     out.int_ci(1,ci) = median(Ym(Yp{ci}>0.9));     % intensity
     out.vol_ci(1,ci) = sum(Yp{ci}(:)>0.5) / 1000;  % volume
  end
  % estimate TIV
  out.TIV = sum(out.vol_ci(1,1:3));

  % create bone mask
  if max(seg8t.lkp) > 6
    Ymsk = smooth3(Yp{4} + Yp{5})>0.5;          % remove small noise
  else
    Ymsk = smooth3(Yp{4})>0.25;                 % remove small noise
  end
  Ymsk = cat_vol_morph(Ymsk,'l');               % get the largest bones element
  Ymsk = cat_vol_morph(Ymsk,'dc',6);            % remove holes in the bone
  
  % extract values (hist)
  out.con       = out.int_ci(1,3);                  % CSF contrast
  out.bone_num  = (sum(Ymsk(:)) / 1000) / out.TIV;  % bone volume normalised by TIV
  out.bone_med  = median(Ym(Ymsk(:))) / out.con;    % median bone intensity value
  out.bone_std  = std(Ym(Ymsk(:)))    / out.con;                     
end
%=======================================================================
function saveSurf(CS,P)
  save(gifti(struct('faces',CS.faces,'vertices',CS.vertices)),P,'Base64Binary'); %#ok<USENS> 
end
%=======================================================================
function CS1 = loadSurf(P)
  CS = gifti(P);
  CS1.vertices = CS.vertices; CS1.faces = CS.faces; 
  if isfield(CS,'cdata'), CS1.cdata = CS.cdata; end
end
%=======================================================================
function x = rgrid(d)
  x = zeros([d(1:3) 3],'single');
  [x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
  for i=1:d(3)
      x(:,:,i,1) = x1;
      x(:,:,i,2) = x2;
      x(:,:,i,3) = single(i);
  end
end
%=======================================================================
function y1 = affind(y0,M)
  y1 = zeros(size(y0),'single');
  for d=1:3
      y1(:,:,:,d) = y0(:,:,:,1)*M(d,1);
      y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,2)*M(d,2);
      y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,3)*M(d,3) + M(d,4);
  end
end
%=======================================================================
function garbage_code
% experement code fragments
% =========================================================================
 if 0
  %% just a table that is nor working
    T = cell2table( num2cell( [ seg8t.mg'; seg8t.mn / tis.seg8o(2); shiftdim(seg8t.vr,4)  ./ tis.seg8o(2)]),...
       'RowName', {'Propotion','Means','Var'}, 'VariableNames', ... 
       cellfun(@num2str,num2cell(seg8t.lkp*1/10 + (1:numel(seg8t.lkp))),'UniformOutput',false) );
     
    uit = uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,... 
      'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0.01 0.90 0.98 0.08]);
  end
 

  if 0
    %%
    Ydiv = cat_vol_div(Ym);
    
    %% get inner and outer bone surfaces 
    Vpp  = cat_io_writenii(Vo, (Ybonepp+1)/3 , '', sprintf('%s.pp','skull') ,...
        'percentage position map - V2 - pial=1/3, central=1/2, white=2/3, 0 and 1 to stablize white/pial values', 'uint8', [0,1/255],...
        min([1 1 2],[1 0 0]),struct());
    cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...   
                   'avg  %0.3f  %0.3f .2  .1  5  0 "0.42"  "0.42"  n 0  0  0 %d %g  0.0 0'], ...          
                    Vpp.fname,out.P.central ,Pouter,-.1, .1, 50, 0.01); 
    cat_system(cmd,opt.verb-3);
    Yppi = max(0,min(1, smooth3(cat_vol_morph(Ybone<.5,'e') + min(0,-Ydiv*4)) + smooth3(Ybonepp))); 
    Vpp  = cat_io_writenii(Vo, Yppi , '', sprintf('%s.pp','skull') ,...
        'percentage position map - V2 - pial=1/3, central=1/2, white=2/3, 0 and 1 to stablize white/pial values', 'uint8', [0,1/255],...
        min([1 1 2],[1 0 0]),struct());
    cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...   
                   'avg  %0.3f  %0.3f .2  .1  5  0 "0.1"  "0.1"  n 0  0  0 %d %g  0.0 0'], ...          
                    Vpp.fname,Pouter,Pouter,-.1, .1, 50, 0.01); 
    cat_system(cmd,opt.verb-3);
    CBSo = loadSurf(Pouter);
  
    %%
    Vpp  = cat_io_writenii(Vo, (Ybonepp+1)/3 , '', sprintf('%s.pp','skull') ,...
        'percentage position map - V2 - pial=1/3, central=1/2, white=2/3, 0 and 1 to stablize white/pial values', 'uint8', [0,1/255],...
        min([1 1 2],[1 0 0]),struct());
    cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...   
                   'avg  %0.3f  %0.3f .2  .1  5  0 "0.57"  "0.57"  n 0  0  0 %d %g  0.0 0'], ...          
                    Vpp.fname,out.P.central ,Pinner,-.1, .1, 50, 0.01); 
    cat_system(cmd,opt.verb-3);
    Yppi = max(0,min(1,Ym + smooth3(Ybrain) + -Ydiv.*(Ybonepp>.5 & ~Ybrain)*10 + (1 - Ybonepp))); 
    Vpp  = cat_io_writenii(Vo, Yppi , '', sprintf('%s.pp','skull') ,...
        'percentage position map - V2 - pial=1/3, central=1/2, white=2/3, 0 and 1 to stablize white/pial values', 'uint8', [0,1/255],...
        min([1 1 2],[1 0 0]),struct());
    cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...   
                   'avg  %0.3f  %0.3f .2  .1  5  0 "0.1"  "0.1"  n 0  0  0 %d %g  0.0 0'], ...          
                    Vpp.fname,Pinner,Pinner,-.1, .1, 50, 0.01); 
    cat_system(cmd,opt.verb-3);
    CBSi = loadSurf(Pinner);
    % get mid-fat surface 
  
    % get outer head surface
  
  %% render
    CBSx = CBSo; maxBM = zeros(size(CBSx.vertices,1),1); minBM = zeros(size(CBSx.vertices,1),1);
    for rpos = 0.0:0.1:1.0
       CBSx = CBSi; CBSx.vertices = CBSi.vertices * rpos  + (1-rpos) * CBSo.vertices; 
       if  abs(rpos - 0.5)<.25
         maxBM = max(maxBM,cat_surf_fun('isocolors',Ym,CBSx, matlab_mm)); 
       end%     else
         minBM = min(minBM,cat_surf_fun('isocolors',Ym,CBSx, matlab_mm)); 
  %     end
    end
    maxBM = maxBM - minBM; 
  
  % S = cat_surf_fun('tlink',CBSi,CBSo);  maxBM = S.facevertexcdata; %maxBM - minBM;S = export(S,'patch');
    %% %Si.vertices = Si.vertices .* repmat(vx_vol,size(Si.vertices,1),1); Si.faces = [Si.faces(:,2) Si.faces(:,1) Si.faces(:,3)];
    CBSx = CBSi; 
    Si.vertices = CBSx.vertices; Si.faces = single(CBSx.faces); Si.facevertexcdata = maxBM; %cat_surf_fun('isocolors',Ym,CBS, matlab_mm); 
    Shi = cat_surf_render2(Si); title(sprintf('bone intensity %s/%s',spm_str_manip(pp,'t'),ff));
    cat_surf_render2('Clim',Shi,[0 1]);  cat_surf_render2('Colorbar',Shi);
  
  end

end
%=======================================================================
