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
  if numel(P)==1
    % if we have only one file the 'C' option of spm_str_manip is not working
    [~,ff]  = spm_fileparts(P{1});
    PC.s    = ff(1);
    PC.m{1} = ff(2:end);
  
    % test for other files
  else
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
    fname    = P{pi}; 
    ind      = max(strfind(spm_fileparts(fname),[filesep 'sub-']));
    der      = max(strfind(spm_fileparts(fname),[filesep 'derivatives' filesep]));
    if ~isempty(der) & ~isempty(ind) 
      derdir{pi}       = fullpath( fileparts(fname(1:ind)) ,'derivatives');  
    else
      derdir{pi}       = ''; 
    end

    sub_ses_anat{pi}   = ''; 
    if ~isempty(ind)
      maindir{pi}      = fileparts(fname(1:ind));  
      sub_ses_anat{pi} = fileparts(fname(ind+1:end));  
    else
      % RD202403:
      % alternative definion based on the depth of the file and is keeping 
      % subdirectories to be more robust in case of a regular but non-BIDS
      % structure wihtout anat directory or with similar filenames, e.g. 
      % for ../derivatives/CAT##.#_#
      %   testdir/subtestdir1/f1.nii
      %   testdir/subtestdir2/f1.nii
      % it result in 
      %   testdir/derivatives/CAT##.#_#/subtestdir1/f1.nii
      %   testdir/derivatives/CAT##.#_#/subtestdir1/f1.nii
      % rather than
      %   testdir/derivatives/CAT##.#_#/f1.nii
      %   testdir/derivatives/CAT##.#_#/f1.nii
      % what would cause conflicts
  
      %%
      subdirs     = strfind(job.output.resdir,['..' filesep]);  
      maindir{pi} = spm_file(fname,'path');
  
      for si = 1:numel(subdirs)
        [maindir{pi},ff,ee] = spm_fileparts(maindir{pi});
        sub_ses_anat{pi} = fullfile([ff ee], sub_ses_anat{pi});
      end
      maindir{pi} = fullfile( maindir{pi},strrep(job.output.resdir,['..' filesep],'')); 
    end
  end

  


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



  for i = 1:numel(P)

    % use subdirectories 
    if job.opts.subdirs
      out(i).P.mridir    = 'mri'; 
      out(i).P.surfdir   = 'surf';
      out(i).P.reportdir = 'report';
    else
      out(i).P.mridir    = ''; 
      out(i).P.surfdir   = ''; 
      out(i).P.reportdir = ''; 
    end
  
    
    % get original file name based on the type of intput segmentation
    if fmethod == 0 || fmethod == 3
      % selection by processed file
      [pp,ff,ee]      = spm_fileparts(P{i}); 
      out(i).CTseg    = contains(ff,'_CTseg');
      if out(i).CTseg
        if strcmp(ff(1:2),'c0'), ffs = 18; elseif strcmp(ff(1:2),'ss'), ffs = 4; end 
        ffe = numel(ff) - 6;
      else
        ffs = prefix1; 
        ffe = numel(ff);
      end
      out(i).P.orgpp  = pp; 
      out(i).P.orgff  = ff(ffs:ffe);
      out(i).P.ppff   = ff;
      out(i).P.ee     = ee; 
      out(i).P.prefix = ff(1:ffs-1);
   
      % input volumes
      out(i).P.org    = fullfile(pp,[ff(ffs:ffe) ee]);
      if out(i).CTseg
        out(i).P.bc   = out(i).P.org;
        out(i).P.seg8 = fullfile(pp,sprintf('mb_fit_CTseg.mat'));
        for ci = 1:5
          out(i).P.cls{ci} = fullfile( pp , sprintf('%s%d%s%s','c0',ci,ff(4:end),ee) ); 
        end
      else % SPM12 case
        out(i).P.bc   = spm_file(out(i).P.org,'prefix','m');
        out(i).P.seg8 = fullfile(pp,sprintf('%s_seg8.mat',out(i).P.orgff));
        for ci = 1:5
          out(i).P.cls{ci} = fullfile( pp , sprintf('%s%d%s%s','c',ci,out(i).P.orgff,ee) ); 
        end
      end
      if ~exist(out(i).P.seg8,'file')
        cat_io_cprintf('err','Cannot process "%s" because the seg8.mat is missing. \n',P{i});
        out(i).process = 0;
      end

    else 
      % selection by RAW file
      [pp,ff,ee]      = spm_fileparts(P{i}); 
      out(i).CTseg    = 0;
      out(i).P.org    = P{i};
      out(i).P.orgpp  = pp; 
      out(i).P.orgff  = ff;
      out(i).P.ee     = ee; 
      out(i).P.ppff   = ['m' ff];
      if job.opts.pmethod == 1
        ffx = 'c'; 
        out(i).P.seg8   = fullfile(pp,sprintf('%s_seg8.mat',out(i).P.orgff));
        out(i).P.bc     = fullfile( pp , sprintf('%s%s%s','m',ff,ee) ); 
        for ci = 1:5
          out(i).P.cls{ci} = fullfile( pp , sprintf('%s%d%s%s',ffx,ci,ff,ee) ); 
        end
      else
        ffx = 'p'; 
        out(i).P.seg8   = fullfile(pp, out(i).P.reportdir, sprintf('cat_%s.xml',out(i).P.orgff));   
        % we have to use the filename from the XML to get the original image 
        % as it could be hidden by BIDS 
        if exist(out(i).P.seg8,'file')
          Sxml          = cat_io_xml( out(i).P.seg8 ); 
          [pp,ff,ee]    = spm_fileparts( Sxml.filedata.fname ); 
          out(i).P.org  = Sxml.filedata.fname;
          out(i).P.bc   = Sxml.filedata.Fm; 
        else
          [pp,ff,ee]    = spm_fileparts(P{i}); 
          out(i).P.org  = P{i};
          out(i).P.bc = fullfile( pp , out(i).P.mridir, sprintf('%s%s%s','m',ff,ee) ); 
        end
        out(i).P.orgpp  = pp; 
        out(i).P.orgff  = ff;
        out(i).P.ee     = ee;
        for ci = 1:5
          out(i).P.cls{ci} = fullfile( pp , out(i).P.mridir, sprintf('%s%d%s%s',ffx,ci,ff,ee) ); 
        end
      end
    end

    % output dirs
    out(i).P.mripath    = fullfile(maindir{pi}, sub_ses_anat{i}, out(i).P.mridir); 
    out(i).P.surfpath   = fullfile(maindir{pi}, sub_ses_anat{i}, out(i).P.surfdir); 
    out(i).P.reportpath = fullfile(maindir{pi}, sub_ses_anat{i}, out(i).P.reportdir); 

    % boney preprocessing mat file for faster reprocessing
    out(i).P.boneyPPmat  = fullfile(out(i).P.mripath, ['boneyPPmat_' ff '.mat']); 
    
    % create dirs if required
    if ~exist(out(i).P.mripath   ,'dir'), mkdir(out(i).P.mripath); end
    if ~exist(out(i).P.surfpath  ,'dir'), mkdir(out(i).P.surfpath); end
    if ~exist(out(i).P.reportpath,'dir'), mkdir(out(i).P.reportpath); end


    % xml/mat output
    out(i).P.report   = fullfile(out(i).P.reportpath, sprintf('%sbonereport%d_%s.jpg', job.output.prefix, job.opts.bmethod, ff));              
    out(i).P.xml      = fullfile(out(i).P.reportpath, sprintf('%s%d_%s.xml'  , job.output.prefix, job.opts.bmethod, ff));
    out(i).P.mat      = fullfile(out(i).P.reportpath, sprintf('%s%d_%s.mat'  , job.output.prefix, job.opts.bmethod, ff));
    
    out(i).P.boneymat = fullfile(out(i).P.reportpath, sprintf('boney%s_%s.mat'  , pmethod, ff));
    out(i).P.boneySPM = fullfile(out(i).P.reportpath, sprintf('boneySPM_%s.mat' , ff));
    out(i).P.boneyCAT = fullfile(out(i).P.reportpath, sprintf('boneyCAT_%s.mat' , ff));

    % vols
    if job.output.writevol
      vols = {'pp','bone','bonethick','headthick'};
      for vi = 1:numel(vols)
        % ... subject affine warped ... 
        prefix = {'','r'}; postfix = {'','_affine'};
        for pi = 1:numel(prefix)
          out(i).P.([prefix{pi} vols{vi} postfix{pi}]) = fullfile(out(i).P.mripath, ...
            sprintf('%s%s%d_%s%s%s', prefix{pi}, vols{vi}, job.opts.bmethod, ff(2:end), postfix{pi}, ee));
        end
      end
    end

    
    % surfs - Are these to be writen anyway?
    if job.opts.bmethod>1 || job.output.writesurf
      out(i).P.central   = fullfile(out(i).P.surfpath, sprintf('%s%d.central.%s.gii', 'bone', job.opts.bmethod, ff));        % central
      out(i).P.marrow    = fullfile(out(i).P.surfpath, sprintf('%s%d.marrow.%s'     , 'bone', job.opts.bmethod, ff));        % marrow
      out(i).P.cortex    = fullfile(out(i).P.surfpath, sprintf('%s%d.cortex.%s'     , 'bone', job.opts.bmethod, ff));        % hard bone - minimal bone values
      out(i).P.thick     = fullfile(out(i).P.surfpath, sprintf('%s%d.thickness.%s'  , 'bone', job.opts.bmethod, ff));        % thickness bone 
      out(i).P.headthick = fullfile(out(i).P.surfpath, sprintf('%s%d.thickness.%s'  , 'head', job.opts.bmethod, ff));        % thickness head
    end
  end
end
