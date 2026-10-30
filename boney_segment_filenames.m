function [out,fmethod,pmethod,P] = boney_segment_filenames(P,job)
%bonefilenames(P). Prepare output filenames and directories. 
%  
% [out,fmethod,pmethod,P] = boney_segment_filenames(P,job)
% 
% P   .. list of raw images or SPM preprocessed images (m*/c*.nii) 
%        (P is updated by removing images that cannot be processed)
% job .. main job structure to include some parameters in the filename etc.
% out .. main output structure
%  .P .. the field with all the filenames
% fmethod .. input type (0 - preprocessed data, 1 - raw data for SPM, 
%            2 - raw data for CAT)
% pmethod .. update preprocessing method 
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________

% TODO: 
% * add CAT preprocessing case (not really relevant as the bone classes are
%   generally not exported) 

  %#ok<*AGROW>
  P = cat_io_strrep(P,',1',''); 

  % == detect the type of input of each file ==
  %  The prefix of preprocessed files is only accepted if the corresponding
  %  SPM/CTseg files exist, to avoid that raw files with similar names 
  %  (e.g. "mprage.nii" or "c1_t1w.nii") are handled as preprocessed data.
  %  ftype .. 0 - raw image, 1 - SPM m*/c* image, 2 - CAT image, 3 - CTseg  
  ftype  = zeros(1,numel(P)); 
  prefix = repmat({''},1,numel(P)); 
  for pi = 1:numel(P)
    [pp,ff] = spm_fileparts(P{pi}); 
    ff = regexprep(ff,'\.nii$',''); % gz
    if contains(ff,'_CTseg') && numel(ff)>2 && any(strcmp(ff(1:2),{'c0','ss'}))
      ftype(pi) = 3; prefix{pi} = ff(1:2); 
    elseif numel(ff)>1 && ff(1)=='m' && ( exist(fullfile(pp,[ff(2:end) '_seg8.mat']),'file') || ...
        exist(fullfile(pp,['c1' ff(2:end) '.nii']),'file') )
      ftype(pi) = 1; prefix{pi} = 'm'; 
    elseif numel(ff)>2 && ~isempty(regexp(ff,'^c[1-6]','once')) && exist(fullfile(pp,[ff(3:end) '_seg8.mat']),'file')
      ftype(pi) = 1; prefix{pi} = ff(1:2); 
    elseif ( numel(ff)>2 && ~isempty(regexp(ff,'^p[0-5]','once')) && exist(fullfile(pp,['m' ff(3:end) '.nii']),'file') ) || ...
           ( numel(ff)>1 && ff(1)=='m' && exist(fullfile(pp,['p1' ff(2:end) '.nii']),'file') )
      ftype(pi) = 2; 
    end
  end
  if any(ftype==2)
    error('boney:segment:CATinput', ['  Preprocessed CAT input (m*/p*) is not supported. \n' ...
      '  Please select the raw images and use the CAT preprocessing option.\n']);
  end
  if numel(unique(ftype))>1
    error('boney:segment:mixedInput', ['  Mixed input of raw images and preprocessed images (m*/c*). \n' ...
      '  Please select only one type.\n']); 
  end

  for pi = 1:numel(P)
  % figure out for each input which is its optimal directory ...
  % - If there is a 'derivatives' directory then we (can) assume it is BIDS 
  %   and that the original files are there.
  % - If there is a 'sub-' directory then we (can) use this to a
  %   derivatives directory.
  % - We also estimate a relative director by the user results directory
  %%
    fname    = P{pi}; 
    subs     = strfind(spm_fileparts(fname),[filesep 'sub-']);
    ders     = strfind(spm_fileparts(fname),[filesep 'derivatives' filesep]);
    sub      = min(subs); % this should not happen
    % in case of multiple derivatives directories we take the one with sub-
    % sibling directories
    der      = ders; 
    for di = 1:numel(ders)
      dersub = cat_vol_findfiles( fileparts(fname(1:ders(di))), 'sub-*', struct('depth',1,'dirs',1,'oneperdir',1));
      if ~isempty(dersub), der = ders(di); end
    end
    subder   = min([sub,der]);

    if isempty(der) && isempty(sub) 
    %% if no sign of BIDS at all then use the user given relative structure
     
      % estimation of user given relative directory 
      maindir = spm_file(fname,'path'); % the orignal 
      maindirr = maindir;               % this is the user specified main home directory for the relative path
      subdirs  = strfind(job.output.resdir,['..' filesep]);  
      sub_ses_anat{pi,1}   = ''; 
      for si = 1:numel(subdirs)
        [maindirr,ff,ee]   = spm_fileparts(maindirr);
        sub_ses_anat{pi,1} = fullfile([ff ee], sub_ses_anat{pi});
      end
      pp_sub_ses_anat{pi,1} = sub_ses_anat{pi,1}; 

      % final relative result directory
      ppdir{pi,1}   = '';
      resdir{pi,1}  = fullfile( maindirr, strrep(job.output.resdir,['..' filesep],'')); 
      sdirs         = max(strfind(job.output.resdir,['..' filesep])); sdirs = sdirs + 1*(sdirs>0); % don't need last seperateor
      outdir{pi,1}  = job.output.resdir(1:sdirs); 
    else 
    %% In this case we have no derivatives dir in the orignal path but we
    %  have the sub- directory that defines it
    
      subinder           = max([subder , subder + min(strfind([filesep fileparts( fname(subder+1:end))],[filesep 'sub-']))]);
      ppdir{pi,1}        = fname(subder+1:subinder-1); 
      resdir{pi,1}       = fullfile( fileparts(fname(1:subder)) , strrep(job.output.resdir,['..' filesep],'')); 
      sub_ses_anat{pi,1} = fileparts(fname(subinder:end)); 
      pp_sub_ses_anat{pi,1} = fileparts(fname(subder+1:end)); 


      %% here we have to find some relative directory discription for the data export 
      if 0 %isempty(der) 
        outdir{pi,1}     = ''; 
      else
        nsdirs           = numel( strfind( pp_sub_ses_anat{pi,1}, filesep) ); nsdirs = nsdirs + (nsdirs>0);
        outdir{pi,1}     = repmat(['..' filesep],1,nsdirs); if numel(outdir{pi})>0, outdir{pi}(end) = []; end % don't need last seperateor
      end
    end

    % tests
    if 0
      %%
      out(pi).P.mridir     = ''; % mri
      out(pi).P.mripath    = fullfile(resdir{pi}, sub_ses_anat{pi}, out(pi).P.mridir); 
      out(pi).P.mrirdir    = fullfile(outdir{pi}, strrep(job.output.resdir,['..' filesep],''), sub_ses_anat{pi}, out(pi).P.mridir); 
      out(pi).P.orgdir     = spm_file(fullfile( fname , out(pi).P.mrirdir ),'fpath'); % removes relative part
      
      fprintf('Dirs %d:\n%s\n',pi,sprintf('  %15s: %s\n', ...
        'boney input',       fname, ...
        'org. dir',          out(pi).P.orgdir, ...
        'boney subdir',      out(pi).P.mripath , ...
        'boney maindir',     resdir{pi}, ...
        'rel sub path',      sub_ses_anat{pi,1}, ...
        'rel pp-sub path',   pp_sub_ses_anat{pi,1}, ...
        'rel pp path',       ppdir{pi})); 
    end
  end
  %%

  if job.opts.subdirs
    mridir    = 'mri'; 
  else
    mridir    = '';
  end

  % basic test for selected preprocessed data m*.nii or c*.nii
  % - in case of preprocessed data the original file is not relevant and
  % - no new preprocessing is done, i.e., missing files have to create errors
  deleteP = false(size(P));
  if ftype(1) == 3
    fmethod = 0; 
    job.opts.pmethod = 3; 
  elseif ftype(1) == 1
    fmethod = 0; 
    if job.opts.pmethod == 2 
      cat_io_cprintf('warn','SPM-preprocessing results are selected although CAT is chosen for processing > Use SPM!\n')
    end
    job.opts.pmethod = 1; 
    for si = 1:numel(P)
      [pp,ff] = spm_fileparts(P{si}); ff = regexprep(ff,'\.nii$','');
      orgff   = ff(numel(prefix{si})+1:end);
      if ~exist( fullfile(pp, sprintf('m%s.nii',orgff)),'file')
        cat_io_cprintf('err',sprintf('  Oops, miss SPM m-file of subject %d.\n',si))
        deleteP(si) = true; continue
      end
      if ~exist( fullfile(pp, sprintf('%s_seg8.mat',orgff)),'file')
        cat_io_cprintf('err',sprintf('  Oops, miss SPM seg8.mat-file of subject %d.\n',si))
        deleteP(si) = true; continue
      end
      for ci = 1:5
        if ~exist( fullfile(pp, sprintf('c%i%s.nii', ci, orgff)),'file')
          cat_io_cprintf('err',sprintf('  Oops, miss SPM-c%d segmentation file of subject %d.\n', ci,si))
          deleteP(si) = true; break
        end
      end
    end
  else
    fmethod = job.opts.pmethod; 
  end 
  P(deleteP) = []; prefix(deleteP) = []; 
  pmethod = job.opts.pmethod; 

  %%
  if isempty(P)
    error(sprintf('  Oops, no files for processing!\n')); 
  end


  for pi = 1:numel(P)

    % use subdirectories 
    if job.opts.subdirs
      out(pi).P.mridir    = 'mri'; 
      out(pi).P.surfdir   = 'surf';
      out(pi).P.reportdir = 'report';
    else
      out(pi).P.mridir    = ''; 
      out(pi).P.surfdir   = ''; 
      out(pi).P.reportdir = ''; 
    end
  
    
    % get original file name based on the type of intput segmentation
    [pp,ff,ee] = spm_fileparts(P{pi}); 
    if strcmp(ee,'.gz') % the derived files are not zipped
      ff = regexprep(ff,'\.nii$',''); ee = '.nii'; 
    end
    if fmethod == 0 
      % selection by processed file
      out(pi).CTseg    = contains(ff,'_CTseg');
      if out(pi).CTseg
        if strcmp(ff(1:2),'c0'), ffs = 18; else, ffs = 4; end 
        ffe = numel(ff) - 6;
      else
        ffs = numel(prefix{pi}) + 1; 
        ffe = numel(ff);
      end
      out(pi).P.orgpp  = pp; 
      out(pi).P.orgff  = ff(ffs:ffe);
      out(pi).P.ppff   = ff;
      out(pi).P.ee     = ee; 
      out(pi).P.prefix = ff(1:ffs-1);
   
      % input volumes
      out(pi).P.org    = fullfile(pp,[ff(ffs:ffe) ee]);
      if ~exist(out(pi).P.org,'file') && exist([out(pi).P.org '.gz'],'file')
        out(pi).P.org  = [out(pi).P.org '.gz']; 
      end
      if out(pi).CTseg
        out(pi).P.bc   = out(pi).P.org;
        out(pi).P.seg8 = fullfile(pp,sprintf('mb_fit_CTseg.mat'));
        for ci = 1:5
          out(pi).P.cls{ci} = fullfile( pp , sprintf('%s%d%s%s','c0',ci,ff(4:end),ee) ); 
        end
      else % SPM12 case
        out(pi).P.bc   = fullfile(pp,['m' out(pi).P.orgff ee]); 
        out(pi).P.seg8 = fullfile(pp,sprintf('%s_seg8.mat',out(pi).P.orgff));
        for ci = 1:numSPMclasses(out(pi).P.seg8,job)
          out(pi).P.cls{ci} = fullfile( pp , sprintf('%s%d%s%s','c',ci,out(pi).P.orgff,ee) ); 
        end
      end

    else 
      % selection by RAW file
      out(pi).CTseg    = 0;
      out(pi).P.org    = P{pi};
      out(pi).P.orgpp  = pp; 
      out(pi).P.orgff  = ff;
      out(pi).P.ee     = ee; 
      out(pi).P.ppff   = ['m' ff];
      if job.opts.pmethod == 1
        ffx = 'c'; 
        out(pi).P.seg8   = fullfile(pp,sprintf('%s_seg8.mat',out(pi).P.orgff));
        out(pi).P.bc     = fullfile( pp , sprintf('%s%s%s','m',ff,ee) ); 
        for ci = 1:numSPMclasses(out(pi).P.seg8,job)
          out(pi).P.cls{ci} = fullfile( pp , sprintf('%s%d%s%s',ffx,ci,ff,ee) ); 
        end
      else
        ffx = 'p'; 
        out(pi).P.seg8   = fullfile(pp, out(pi).P.reportdir, sprintf('cat_%s.xml',out(pi).P.orgff));   
        % we have to use the filename from the XML to get the original image 
        % as it could be hidden by BIDS 
        if exist(out(pi).P.seg8,'file')
          Sxml          = cat_io_xml( out(pi).P.seg8 ); 
          [pp,ff,ee]    = spm_fileparts( Sxml.filedata.fname ); 
          out(pi).P.org  = Sxml.filedata.fname;
          out(pi).P.bc   = Sxml.filedata.Fm; 
        else
          out(pi).P.org  = P{pi};
          out(pi).P.bc   = fullfile( pp , out(pi).P.mridir, sprintf('%s%s%s','m',ff,ee) ); 
        end
        out(pi).P.y      = fullfile( pp , out(pi).P.mridir, sprintf('%s%s%s','y_',ff,ee) ); 
        out(pi).P.orgpp  = pp; 
        out(pi).P.orgff  = ff;
        out(pi).P.ee     = ee;
        for ci = 1:5
          out(pi).P.cls{ci} = fullfile( pp , out(pi).P.mridir, sprintf('%s%d%s%s',ffx,ci,ff,ee) ); 
        end
      end
    end

    
    %% output dirs
    out(pi).P.resdir     = resdir{pi};
    out(pi).P.mripath    = fullfile(resdir{pi}, sub_ses_anat{pi}, out(pi).P.mridir); 
    out(pi).P.surfpath   = fullfile(resdir{pi}, sub_ses_anat{pi}, out(pi).P.surfdir); 
    out(pi).P.reportpath = fullfile(resdir{pi}, sub_ses_anat{pi}, out(pi).P.reportdir); 

    % output dirs for CAT volume export fuction 
    out(pi).P.mrirdir    = fullfile(outdir{pi}, strrep(job.output.resdir,['..' filesep],''), sub_ses_anat{pi}, out(pi).P.mridir); 
    out(pi).P.surfrdir   = fullfile(outdir{pi}, strrep(job.output.resdir,['..' filesep],''), sub_ses_anat{pi}, out(pi).P.surfdir); 
    out(pi).P.reportrdir = fullfile(outdir{pi}, strrep(job.output.resdir,['..' filesep],''), sub_ses_anat{pi}, out(pi).P.reportdir); 
    out(pi).P.orgdir     = spm_file(fullfile( out(pi).P.orgpp , out(pi).P.mrirdir),'fpath'); 
    
    
    %% boney preprocessing mat file for faster reprocessing
    PP = {'SPM12','CAT12','CTseg'}; 
    out(pi).P.boneyPPmat  = fullfile(out(pi).P.mripath, ['boney' PP{job.opts.pmethod} 'PPmat_' out(pi).P.orgff '.mat']); 
    
    % create dirs if required
    if ~exist(out(pi).P.mripath   ,'dir'), mkdir(out(pi).P.mripath); end
    if ~exist(out(pi).P.surfpath  ,'dir'), mkdir(out(pi).P.surfpath); end
    if ~exist(out(pi).P.reportpath,'dir'), mkdir(out(pi).P.reportpath); end


    % xml/mat output
    out(pi).P.report   = fullfile(out(pi).P.reportpath, sprintf('%sbonereport_%s.jpg', job.output.prefix, out(pi).P.orgff));              
    out(pi).P.xml      = fullfile(out(pi).P.reportpath, sprintf('%s%s.xml'          , job.output.prefix, out(pi).P.orgff));
    out(pi).P.mat      = fullfile(out(pi).P.reportpath, sprintf('%s%s.mat'          , job.output.prefix, out(pi).P.orgff));
    
    out(pi).P.boneymat = fullfile(out(pi).P.reportpath, sprintf('boney_%s.mat'    , out(pi).P.orgff));
    out(pi).P.boneySPM = fullfile(out(pi).P.reportpath, sprintf('boneySPM_%s.mat' , out(pi).P.orgff));
    out(pi).P.boneyCAT = fullfile(out(pi).P.reportpath, sprintf('boneyCAT_%s.mat' , out(pi).P.orgff));
    out(pi).P.tpmsel   = fullfile(out(pi).P.reportpath, sprintf('boneyTPM_%s.mat' , out(pi).P.orgff)); % automatic TPM selection

    % vols (names as written by cat_io_writenii in boney_segment_extractbone)
    if job.output.writevol
      vols = {'bonepp','bone','bonethick','headthick'};
      vpre = {'bonepp','bonemarrow','bonethickness','headthickness'};
      for vi = 1:numel(vols)
        out(pi).P.(vols{vi})                  = fullfile(out(pi).P.mripath, sprintf('%s_%s%s'  , vpre{vi}, out(pi).P.orgff, ee));
        out(pi).P.(['r' vols{vi} '_affine'])  = fullfile(out(pi).P.mripath, sprintf('r%s_%s_affine%s', vpre{vi}, out(pi).P.orgff, ee));
      end
    end

    
    % surfs - Are these to be writen anyway?
    if job.opts.bmethod>1 || job.output.writesurf
      out(pi).P.central   = fullfile(out(pi).P.surfpath, sprintf('%s.central.%s.gii', 'bone', out(pi).P.orgff));        % central
      out(pi).P.marrow    = fullfile(out(pi).P.surfpath, sprintf('%s.marrow.%s'     , 'bone', out(pi).P.orgff));        % marrow
      out(pi).P.cortex    = fullfile(out(pi).P.surfpath, sprintf('%s.cortex.%s'     , 'bone', out(pi).P.orgff));        % hard bone - minimal bone values
      out(pi).P.thick     = fullfile(out(pi).P.surfpath, sprintf('%s.thickness.%s'  , 'bone', out(pi).P.orgff));        % thickness bone 
      out(pi).P.headthick = fullfile(out(pi).P.surfpath, sprintf('%s.thickness.%s'  , 'head', out(pi).P.orgff));        % thickness head
    end
  end

end
%=======================================================================
function ncls = numSPMclasses(Pseg8,job)
%numSPMclasses. Number of written SPM classes (6 for 7-class TPMs such as
% the BlaiottaTPM) that is defined by the existing seg8 file or by the
% selected TPM.
  ncls = 5;
  if exist(Pseg8,'file')
    try
      s8   = load(Pseg8,'lkp');
      ncls = 5 + (max(s8.lkp)==7);
    end
  else
    Ptpms = boney_segment_preprocessing('tpms');
    if ~isempty(job.opts.spm.tpm)
      Ptpm = spm_file(job.opts.spm.tpm,'number','');
    elseif job.opts.ctpm > 0
      Ptpm = Ptpms{job.opts.ctpm};
    else
      Ptpm = '';
    end
    if ~isempty(Ptpm) && exist(Ptpm,'file') && numel(spm_vol(Ptpm)) == 7
      ncls = 6;
    end
  end
end
