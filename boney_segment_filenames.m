function out = boney_segment_filenames(P,job)
%bonefilenames(P). Prepare output filesnames and directories. 
%  
% out = boney_segment_filenames(P,job)
% 
% P   .. list of m-files
% job .. main job structure to include some parameters in the filename etc.
% out .. main output structure
%  .P .. the field with all the filenames
% _________________________________________________________________________
%
% Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________
% $Id$

% TODO: 
% * add CAT preprocessing case
% * creation of surface names allways required?

  %#ok<*AGROW> 

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
    else % SPM case
      out(i).P.bc   = P{i};
      out(i).P.seg8 = fullfile(pp,sprintf('%s_seg8.mat',out(i).P.orgff));
% ############## possible CAT case that would have to use the CAT xml/mat file      
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