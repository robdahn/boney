function [pp,job,stopit] = boney_segment_checkPP(out,job)
%boney_segment_checkPP. Check existing SPM/CAT preprocessing results.
% Test for each subject if the SPM/CAT segmentation is available and if
% SPM segmentations were processed with the current SPM parameters (TPM,
% number of Gaussians, bias FWHM, sampling distance). Reusing data is
% always risky, so differences are reported and the user decides how to
% proceed (only in the main process and if the GUI is available):
%   - use existing segmentation (default)
%   - rerun the segmentation with the new SPM parameters
%   - stop processing (e.g. to use another result directory)
% Without GUI, the opts.prerun setting is used (see boney_segment).
%
%  [pp,job,stopit] = boney_segment_checkPP(out,job)
%
%  out        .. output structure with filenames (see boney_segment_filenames)
%  job        .. main job structure (opts.prerun is updated by the decision
%                of the user to avoid further questions)
%  pp         .. structure with logical vectors
%   .proc     .. (re)run segmentation (missing data or user/parameter choice)
%   .force    .. existing results have to be updated (no lazy processing)
%   .bconly   .. only the bias corrected image m* is missing
%   .diff     .. cell with the parameter differences of each subject
%  stopit     .. stop processing
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________

  n         = numel(out);
  pp.proc   = false(1,n);
  pp.force  = false(1,n);
  pp.bconly = false(1,n);
  pp.diff   = cell(1,n);
  stopit    = false;

  % processed input (m*/c*) - nothing to do
  if ~job.opts.fmethod, return; end

  for i = 1:n
    if job.opts.pmethod == 2
      % CAT: only test if the files exist (the CAT parameters are not comparable)
      pp.proc(i) = ~exist(out(i).P.cls{1},'file') || ~exist(out(i).P.bc,'file') || ~exist(out(i).P.seg8,'file');
    else
      % SPM: test files and parameters
      if ~exist(out(i).P.seg8,'file') || ~all(cellfun(@(f) exist(f,'file')>0, out(i).P.cls))
        pp.proc(i) = true;
      else
        pp.bconly(i) = ~exist(out(i).P.bc,'file');
        pp.diff{i}   = compareSPMparameter(out(i).P.seg8, job);
      end
    end
  end

  % always rerun preprocessing
  if job.opts.prerun == 1
    pp.proc(:)   = true;
    pp.force(:)  = true;
    pp.bconly(:) = false;
    return
  end


  % handle differences in SPM parameters
  ndiff = ~cellfun('isempty',pp.diff);
  if any(ndiff)
    % print differences
    if job.opts.verb
      cat_io_cprintf('warn',sprintf(['\n  The SPM parameters of %d of %d existing segmentations differ from the current setting:\n'], ...
        sum(ndiff), n));
      fi = find(ndiff);
      for i = fi(1:min(numel(fi),10))
        cat_io_cprintf('warn',sprintf('    %s: %s\n', spm_str_manip(out(i).P.org,'a50'), strjoin(pp.diff{i},', ')));
      end
      if numel(fi)>10, cat_io_cprintf('warn',sprintf('    ... and %d further cases\n', numel(fi) - 10)); end
    end

    % decision
    prerun = job.opts.prerun;
    if prerun == 0
      if ~isfield(job,'process_index') && usejava('desktop') && feature('ShowFigureWindows') && ~isempty(spm_figure('FindWin','Interactive'))
        sel = spm_input(sprintf('SPM parameters differ in %d segmentation(s):',sum(ndiff)), 1, 'm', ...
          'use existing segmentation|rerun segmentation with new SPM parameters|stop processing', [2 3 0], 1);
        if sel == 0
          stopit = true;
          if job.opts.verb, cat_io_cprintf('warn','  Stopped by user (e.g. to use another result directory).\n\n'); end
          return
        end
        prerun = sel;
      else
        prerun = 2; % use existing segmentation
      end
    end
    if prerun == 3
      pp.proc(ndiff)   = true;
      pp.force(ndiff)  = true;
      pp.bconly(ndiff) = false;
      if job.opts.verb, cat_io_cprintf('warn','  Rerun segmentation with the new SPM parameters.\n\n'); end
    else
      if job.opts.verb, cat_io_cprintf('warn','  Use existing segmentation(s).\n\n'); end
    end
    job.opts.prerun = prerun; % avoid further questions (e.g. in parallel processing)
  end
end
%=======================================================================
function sdiff = compareSPMparameter(Pseg8, job)
%compareSPMparameter. Compare the parameters of an existing SPM
% segmentation with the current settings. The bias FWHM and the sampling
% distance are not stored in the seg8 file but are defined by the size of
% the bias field (Tbias) and the subsampling matrix (MT) (see spm_preproc8).

  sdiff = {};
  try
    s8 = load(Pseg8,'tpm','lkp','image','Tbias','MT');
  catch
    sdiff = {'cannot read seg8.mat'};
    return
  end

  % TPM
  [~,tpmff] = spm_fileparts(s8.tpm(1).fname);
  Ptpms     = boney_segment_preprocessing('tpms');
  if ~isempty(job.opts.spm.tpm)
    [~,exptpm] = spm_fileparts(job.opts.spm.tpm); exptpm = {exptpm};
  elseif job.opts.ctpm == 0 % automatic selection and fallback, i.e., all TPMs are fine
    exptpm = cellfun(@(f) spm_str_manip(f,'tr'), Ptpms, 'UniformOutput', false);
  else
    exptpm = {spm_str_manip(Ptpms{job.opts.ctpm},'tr')};
  end
  if ~any(strcmp(tpmff,exptpm))
    sdiff{end+1} = sprintf('TPM %s (not %s)',tpmff,strjoin(exptpm,'/'));
  end

  % Gaussians
  ngaus = histcounts(s8.lkp,0.5:1:max(s8.lkp)+0.5);
  if isempty(job.opts.spm.ngaus)
    expngaus = boney_segment_preprocessing('ngaus',s8.tpm(1).fname);
  else
    expngaus = job.opts.spm.ngaus;
  end
  if numel(ngaus) ~= numel(expngaus) || any(ngaus ~= expngaus)
    sdiff{end+1} = sprintf('ngaus [%s] (not [%s])',num2str(ngaus),num2str(expngaus));
  end

  % bias FWHM: d3 = ceil(2 * vx * d0 / fwhm)
  vx = sqrt(sum(s8.image(1).mat(1:3,1:3).^2));
  d0 = s8.image(1).dim(1:3);
  if isempty(job.opts.spm.biasfwhm)
    expfwhm = 60 - 20*job.opts.bias;
  else
    expfwhm = job.opts.spm.biasfwhm;
  end
  if isfield(s8,'Tbias') && ~isempty(s8.Tbias) && ~isinf(expfwhm)
    d3 = size(s8.Tbias{1}); d3(end+1:3) = 1;
    if any( d3 ~= ceil(2 * vx .* d0 / expfwhm) )
      sdiff{end+1} = sprintf('biasfwhm ~%0.0f mm (not %0.0f mm)', mean(2 * vx .* d0 ./ d3), expfwhm);
    end
  end

  % sampling distance: sk = round(samp/vx)
  if isfield(s8,'MT') && ~isempty(s8.MT)
    sk = diag(s8.MT(1:3,1:3))';
    if any( sk ~= max([1 1 1],round(job.opts.spm.samp ./ vx)) )
      sdiff{end+1} = sprintf('samp ~%0.1f mm (not %0.1f mm)', mean(sk .* vx), job.opts.spm.samp);
    end
  end
end
