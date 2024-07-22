function [out,fmethod,pmethod] = boney_segment_filenames(P,job)
%bonefilenames(P). Prepare output filenames and directories. 
%  
% [out,fmethod,pmethod] = boney_segment_filenames(P,job)
% 
% P   .. list of m-files
% job .. main job structure to include some parameters in the filename etc.
% out .. main output structure
%  .P .. the field with all the filenames
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

  % try to get the prefix
  P = cat_io_strrep(P,',1',''); 
  if isscalar(P)
    % if we have only one file the 'C' option of spm_str_manip is not working
    [~,ff]  = spm_fileparts(P{1});
    PC.s    = ff(1);
    PC.m{1} = ff(2:end);
  else
    % test for other files
    [~,PC]  = spm_str_manip(P,'tC'); 
    if numel(PC.s)>1 && PC.s(1)=='m'
      for mi=1:numel(PC.m)
        PC.m{mi} = [PC.s(2:end) PC.m{mi}];
      end
      PC.s = PC.s(1);
    end
    if contains(PC.e,'CTseg')
      for mi=1:numel(PC.m)
        PC.m{mi} = [PC.s(4:end) PC.m{mi}];
      end
      PC.s = PC.s(1:3);
    end
    if numel(PC.e)>4 
      for mi=1:numel(PC.m)
        PC.m{mi} = [PC.m{mi} PC.e(1:end-4)];
      end
      PC.e = PC.e(end-3:end);
    end
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
      dersub = cat_vol_findfiles( fileparts(fname(1:ders)), 'sub-*', struct('depth',1,'dirs',1,'oneperdir',1));
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
    if 1
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

  % basic test for selected preprocessed data m*.nii, c*.nii, or p*.nii 
  % - in case of preprocessed data the original file is not relevant and
  % - no new preprocessing is done, i.e., missing files have to create errors
  % - CAT supports BIDS that makes it complicated
  fmethod = 0; out = [];  
  if ~isempty(PC.s) && ( strcmp(PC.s(1),'m') || strcmp(PC.s(1),'c') || strcmp(PC.s(1),'p') ) 
    if strcmp(PC.s(1),'m')
      % try to handle the m-file input by updating the pmethod 
      prefix1 = 2; 
      if job.opts.pmethod == 0 && ...
        exist( fullfile(fileparts(P{1}), ['c1' PC.m{1} '.nii']),'file') && ...
        exist( fullfile(fileparts(P{1}), mridir, ['p1' PC.m{1} '.nii']),'file') 
        % if m* input is used and SPM and CAT segments are available and 
        % neither SPM nor CAT is selected, the user has to select a
        % segmentation
          error('boney:segment:unclearSelection', ...
            '  Oops, should I use SPM or CAT? Please select the SPM-c1 or the CAT-p1 segmentation files!\n')
      else
        % if m* input is used and we have either SPM or CAT and previous 
        % segmentation flag is used, then update based on the input
        % segmentation 
        for si = 1:numel(P)
          for ci = 1:5
            if job.opts.pmethod == 0
              if exist( fullfile(fileparts(P{si}), sprintf('%s%i%s.nii','c', ci, PC.m{1})),'file') 
                job.opts.pmethod = 1; 
              elseif exist( fullfile(fileparts(P{si}), sprintf('%s%i%s.nii','p1', ci, PC.m{1})),'file') 
                job.opts.pmethod = 2;
              else
                error('boney:segment:missingSegmentationInput', ...
                  '  Oops, miss SPM-c%d and/or CAT-p%d segmentation file of subject %d.\n',ci,ci,si);
              end
            else
              if exist( fullfile(fileparts(P{si}), sprintf('%s%i%s.nii','c', ci, PC.m{1})),'file') && job.opts.pmethod==2
                error('boney:segment:unclearSetup', ...
                  '  Oops, found SPM c*-segmentation files but you selected CAT segmentation.\n');
              elseif exist( fullfile(fileparts(P{si}), sprintf('%s%i%s.nii','p1', ci, PC.m{1})),'file') && job.opts.pmethod==1
                error('boney:segment:unclearSetup', ...
                  '  Oops, found CAT p*-segmentation files but you selected SPM segmentation.\n');
              end
            end
          end
        end
      end
    elseif strcmp(PC.s(1:2),'c0')  % CTseg
      prefix1 = 3; 
      fmethod = 3;
      job.opts.pmethod = 3; 
      if job.opts.pmethod == 1 
        cat_io_cprintf('warn','CTseg-preprocessing results are selected although SPM is chosen for processing > Use CTseg!\n')
      elseif job.opts.pmethod == 2 
        cat_io_cprintf('warn','CTseg-preprocessing results are selected although CAT is chosen for processing > Use CTseg!\n')
      end
      job.opts.pmethod = 3; 
    else
      prefix1 = 3; 
      % if SPM or CAT segments are selected then just update the setting
      fmethod = 1 + PC.s(1)=='p'; 
      if all([ fmethod job.opts.pmethod ] == [ 1 2 ])
        cat_io_cprintf('warn','SPM-preprocessing results are selected although CAT is chosen for processing > Use SPM!\n')
        job.opts.pmethod = 1; 
      elseif all([ fmethod job.opts.pmethod ] == [ 2 1 ])
        cat_io_cprintf('warn','CAT-preprocessing results are selected although SPM is chosen for processing > Use CAT!\n')
        job.opts.pmethod = 2; 
      end
    end

    %% test for preprocessed data
    deleteP = false(size(P));
    switch job.opts.pmethod 
      case 1, pmethod = 'SPM';   sprefix = 'c';  
      case 2, pmethod = 'CAT';   sprefix = 'p';  
      case 3, pmethod = 'CTseg'; sprefix = 'c0'; 
    end
    for si = 1:numel(P)
      if job.opts.pmethod~=3 && ~exist( fullfile(fileparts(P{si}), sprintf('m%s.nii',PC.m{si})),'file')
        cat_io_cprintf('err',sprintf('  Oops, miss SPM/CAT m-file of subject %d.\n',si))
        deleteP(si) = true; 
      end
      for ci = 1:5
        if ~exist( fullfile(fileparts(P{si}), sprintf('%s%i%s.nii',sprefix,  ci, PC.m{si})),'file')
          cat_io_cprintf('err',sprintf('  Oops, miss %s-%s%d segmentation file of subject %d.\n', ...
            pmethod,sprefix,ci,si))
          deleteP(si) = true; 
        end
      end
    end
    P(deleteP) = [];
  else
    fmethod = job.opts.pmethod; 
    prefix1 = 1; 
  end 
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
    if fmethod == 0 || fmethod == 3
      % selection by processed file
      [pp,ff,ee]      = spm_fileparts(P{pi}); 
      out(pi).CTseg    = contains(ff,'_CTseg');
      if out(pi).CTseg
        if strcmp(ff(1:2),'c0'), ffs = 18; elseif strcmp(ff(1:2),'ss'), ffs = 4; end 
        ffe = numel(ff) - 6;
      else
        ffs = prefix1; 
        ffe = numel(ff);
      end
      out(pi).P.orgpp  = pp; 
      out(pi).P.orgff  = ff(ffs:ffe);
      out(pi).P.ppff   = ff;
      out(pi).P.ee     = ee; 
      out(pi).P.prefix = ff(1:ffs-1);
   
      % input volumes
      out(pi).P.org    = fullfile(pp,[ff(ffs:ffe) ee]);
      if out(pi).CTseg
        out(pi).P.bc   = out(pi).P.org;
        out(pi).P.seg8 = fullfile(pp,sprintf('mb_fit_CTseg.mat'));
        for ci = 1:5
          out(pi).P.cls{ci} = fullfile( pp , sprintf('%s%d%s%s','c0',ci,ff(4:end),ee) ); 
        end
      else % SPM12 case
        out(pi).P.bc   = spm_file(out(pi).P.org,'prefix','m');
        out(pi).P.seg8 = fullfile(pp,sprintf('%s_seg8.mat',out(pi).P.orgff));
        for ci = 1:5
          out(pi).P.cls{ci} = fullfile( pp , sprintf('%s%d%s%s','c',ci,out(pi).P.orgff,ee) ); 
        end
      end
      if ~exist(out(pi).P.seg8,'file')
        cat_io_cprintf('err','Cannot process "%s" because the seg8.mat is missing. \n',P{pi});
        out(pi).process = 0;
      end

    else 
      % selection by RAW file
      [pp,ff,ee]      = spm_fileparts(P{pi}); 
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
        for ci = 1:5
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
          [pp,ff,ee]    = spm_fileparts(P{pi}); 
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
    out(pi).P.boneyPPmat  = fullfile(out(pi).P.mripath, ['boney' PP{job.opts.pmethod} 'PPmat_' ff '.mat']); 
    
    % create dirs if required
    if ~exist(out(pi).P.mripath   ,'dir'), mkdir(out(pi).P.mripath); end
    if ~exist(out(pi).P.surfpath  ,'dir'), mkdir(out(pi).P.surfpath); end
    if ~exist(out(pi).P.reportpath,'dir'), mkdir(out(pi).P.reportpath); end


    % xml/mat output
    out(pi).P.report   = fullfile(out(pi).P.reportpath, sprintf('%sbonereport_%s.jpg', job.output.prefix, ff));              
    out(pi).P.xml      = fullfile(out(pi).P.reportpath, sprintf('%s%s.xml'          , job.output.prefix, ff));
    out(pi).P.mat      = fullfile(out(pi).P.reportpath, sprintf('%s%s.mat'          , job.output.prefix, ff));
    
    out(pi).P.boneymat = fullfile(out(pi).P.reportpath, sprintf('boney_%s.mat'    , ff));
    out(pi).P.boneySPM = fullfile(out(pi).P.reportpath, sprintf('boneySPM_%s.mat' , ff));
    out(pi).P.boneyCAT = fullfile(out(pi).P.reportpath, sprintf('boneyCAT_%s.mat' , ff));

    % vols
    if job.output.writevol
      vols = {'pp','bone','bonethick','headthick'};
      for vi = 1:numel(vols)
        % ... subject affine warped ... 
        prefix = {'','r'}; postfix = {'','_affine'};
        for pri = 1:numel(prefix)
          out(pi).P.([prefix{pri} vols{vi} postfix{pri}]) = fullfile(out(pi).P.mripath, ...
            sprintf('%s%s_%s%s%s', prefix{pri}, vols{vi}, ff(2:end), postfix{pri}, ee));
        end
      end
    end

    
    % surfs - Are these to be writen anyway?
    if job.opts.bmethod>1 || job.output.writesurf
      out(pi).P.central   = fullfile(out(pi).P.surfpath, sprintf('%s.central.%s.gii', 'bone', ff));        % central
      out(pi).P.marrow    = fullfile(out(pi).P.surfpath, sprintf('%s.marrow.%s'     , 'bone', ff));        % marrow
      out(pi).P.cortex    = fullfile(out(pi).P.surfpath, sprintf('%s.cortex.%s'     , 'bone', ff));        % hard bone - minimal bone values
      out(pi).P.thick     = fullfile(out(pi).P.surfpath, sprintf('%s.thickness.%s'  , 'bone', ff));        % thickness bone 
      out(pi).P.headthick = fullfile(out(pi).P.surfpath, sprintf('%s.thickness.%s'  , 'head', ff));        % thickness head
    end
  end

  if 0
    out(pi).P
  end
end
