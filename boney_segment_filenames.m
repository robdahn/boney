function [out,fmethod] = boney_segment_filenames(P,job)
%bonefilenames(P). Prepare output filesnames and directories. 
%  
% [out,pmethod] = boney_segment_filenames(P,job)
% 
% P   .. list of m-files
% job .. main job structure to include some parameters in the filename etc.
% out .. main output structure
%  .P .. the field with all the filenames
% pmethod .. update preprocessing method 
% _________________________________________________________________________
%
% Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________

% TODO: 
% * add CAT preprocessing case
% * creation of surface names allways required?

  %#ok<*AGROW>

  % try to get the prefix
  if numel(P)==1
    % if we have only one file the 'C' option of spm_str_manip is not working
    [~,ff]  = spm_fileparts(P{1});
    PC.s    = ff(1);
    PC.m{1} = ff(2:end);
  
    % test for other files
  else
    [~,PC]  = spm_str_manip(P,'tC'); 
  end

  % basic test for selected preprocessed data m*.nii, c*.nii, or p*.nii 
  % - in case of preprocessed data the original file is not relevant and
  % - no new preprocessing is done, i.e., missing files have to create errors
  % - CAT support BIDS what makes it complicated
  fmethod = 0; 
  if ~isempty(PC.s) && ( strcmp(PC.s(1),'m') || strcmp(PC.s(1),'c') || strcmp(PC.s(1),'p') ) 
    if strcmp(PC.s(1),'m')
      % try to handle the m-file input by updating the pmethod 
      prefix1 = 2; 
      if job.opts.pmethod == 0 && ...
        exist( fullfile(fileparts(P{1}), ['c1' PC.m{1} '.nii']),'file') && ...
        exist( fullfile(fileparts(P{1}), ['p1' PC.m{1} '.nii']),'file') 
        % if m* input is used and SPM and CAT segments are available and 
        % neither SPM nor CAT is selected the uses has to select a
        % segmentation
          error('Ups, should I use SPM or CAT? Please select the SPM-c1 or the CAT-p1 segmentation files!')
      else
        % if m* input is used and we have either SPM or CAT and previous 
        % segmentation flag is used than update based on the input
        % segmentation 
        for si = 1:numel(P)
          for ci = 1:5
            if exist( fullfile(fileparts(P{si}), sprintf('%s%i%s.nii','c', ci, PC.m{1})),'file')
              job.opts.pmethod = 1; fmethod = 1; 
            elseif exist( fullfile(fileparts(P{si}), sprintf('%s%i%s.nii','p1', ci, PC.m{1})),'file')
              job.opts.pmethod = 2; fmethod = 2; 
            else
              error('Ups, miss SPM-c1 and/or CAT-p%d segmentation file of subject %d: ',ci,si)
            end
          end
        end
      end
    else
      prefix1 = 3; 
      % if SPM or CAT segments are selected then just update the setting
      fmethod = 1 + PC.s(1)=='p'; 
      switch [ fmethod job.opts.pmethod ]
        case [ 1 2 ]
          cat_io_cprintf('warn','SPM-preprocessing results are selected alhtough CAT is choosen for processing > Use SPM!\n')
          job.opts.pmethod = 1; 
        case [ 2 1 ]
          cat_io_cprintf('warn','CAT-preprocessing results are selected alhtough SPM is choosen for processing > Use CAT!\n')
          job.opts.pmethod = 2; 
      end
    end

    % test for preprocessed data
    deleteP = false(size(P)); 
    for si = 1:numel(P)
      for ci = 1:5
        if ~exist( fullfile(fileparts(P{si}), sprintf('%s%i%s.nii',PC.s(1),  ci, PC.m{1})),'file')
          cat_io_cprintf('err','Ups, miss SPM-c1 and/or CAT-p%d segmentation file of subject %d: ',ci,si)
          deleteP(si) = true; 
        end
      end
    end
    P(deleteP) = [];
  else
    prefix1 = 1; 
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
    if fmethod == 0 
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
      else % SPM case
        out(i).P.bc   = P{i};
        out(i).P.seg8 = fullfile(pp,sprintf('%s_seg8.mat',out(i).P.orgff));
  % ############## possible CAT case that would have to use the CAT xml/mat file      
      end
      if ~exist(out(i).P.seg8,'file')
        cat_io_cprintf('err','Cannot process "%s" because the seg8.mat is missing. \n',P{i});
        out(i).process = 0;
      end

    else 
      [pp,ff,ee]      = spm_fileparts(P{i}); 
      out(i).CTseg    = 0;
      out(i).P.org    = P{i};
      out(i).P.orgpp  = pp; 
      out(i).P.orgff  = ff;
      out(i).P.ee     = ee; 
      out(i).P.bc     = fullfile( pp , sprintf('%s%s%s','m',ff,ee) ); 
      out(i).P.ppff   = ['m' ff];
      %out(i).P.p0     = ['p0' ff];
      if job.opts.pmethod==1
        ffx = 'c'; 
        out(i).P.seg8 = fullfile(pp,sprintf('%s_seg8.mat',out(i).P.orgff));
      else
        ffx = 'p'; 
        out(i).P.seg8 = fullfile(pp,sprintf('cat_%s.xml',out(i).P.orgff));   
      end
      for ci = 1:5
        out(i).P.cls{ci} = fullfile( pp , sprintf('%s%d%s%s',ffx,ci,ff,ee) ); 
      end
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
    out(i).P.report = fullfile(out(i).P.reportpath, sprintf('bonereport%d_%s.jpg', job.opts.bmethod, ff(2:end)));              
    out(i).P.xml    = fullfile(out(i).P.reportpath, sprintf('catbones%d_%s.xml'  , job.opts.bmethod, ff(2:end)));
    out(i).P.mat    = fullfile(out(i).P.reportpath, sprintf('catbones%d_%s.mat'  , job.opts.bmethod, ff(2:end)));
    

    % vols
    if job.output.writevol
      vols = {'pp','bonemarrow','headthick'};
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
% #################    
    if 1 %job.output.writesurf
      out(i).P.central   = fullfile(out(i).P.surfpath, sprintf('%s%d.central.%s.gii', 'bone', job.opts.bmethod, ff));        % central
      out(i).P.marrow    = fullfile(out(i).P.surfpath, sprintf('%s%d.marrow.%s'     , 'bone', job.opts.bmethod, ff));        % marrow
      out(i).P.thick     = fullfile(out(i).P.surfpath, sprintf('%s%d.thickness.%s'  , 'bone', job.opts.bmethod, ff));        % thickness
      out(i).P.headthick = fullfile(out(i).P.surfpath, sprintf('%s%d.thickness.%s'  , 'head', job.opts.bmethod, ff));        % thickness
      out(i).P.marrowmin = fullfile(out(i).P.surfpath, sprintf('%s%d.marrowmin.%s'  , 'bone', job.opts.bmethod, ff));        % minimal bone values
      out(i).P.marrowmax = fullfile(out(i).P.surfpath, sprintf('%s%d.marrowmax.%s'  , 'bone', job.opts.bmethod, ff));        % maximum bone values
    end
  end
end