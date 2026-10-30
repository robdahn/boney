function [ctpm, sel] = boney_segment_selectTPM(Pimg, verb)
%boney_segment_selectTPM. Automatic selection of the TPM for SPM.
% The standard SPM TPM is not optimal for children (very thin skull, small
% head) and the elderly (thicker and fatty skull, larger head), as the head
% and skull change strongly with age. Each candidate TPM is registered by a
% fast affine registration (spm_maff8 with low sampling) to estimate how
% good the TPM fits to the image (mutual information). In addition, the
% scaling of the affine registration (head size relative to the TPM) is
% used. Both are combined into an age-like parameter that is only used for
% the selection (and printed for humans) but that is not a measure.
%
%  [ctpm, sel] = boney_segment_selectTPM(Pimg, verb)
%
%  Pimg       .. image filename (not zipped)
%  verb       .. print selection (default 0)
%  ctpm       .. selected TPM (1-SPM TPM, 2-children TPM, 3-BlaiottaTPM)
%  sel        .. structure with selection details
%   .Ptpm     .. tested TPMs
%   .ll       .. log-likelihood (mutual information) of the affine
%                 registration of each TPM
%   .dll      .. difference of the log-likelihood to the SPM TPM
%   .zoom     .. isotropic scaling of the affine registration to the SPM
%                 TPM (>1 for smaller heads)
%   .agelike  .. age-like parameter (only for humans)
%   .ctpm     .. selected TPM
%   .time     .. processing time
%   .fallback .. originally selected TPM if the segmentation failed and
%                was repeated with another TPM (0 - no fallback; set by
%                boney_segment)
%   .tried    .. all used TPMs (selection and fallbacks)
%
% See also boney_segment_preprocessing.
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________

  if ~exist('verb','var'), verb = 0; end
  stime = clock;

  % thresholds of the selection (boney development set, 66 scans):
  %  The SPM TPM gave the best test-retest and protocol consistency and
  %  the best fat measures, so the other TPMs have to fit clearly better.
  %  The BlaiottaTPM (extra fat and soft tissue classes) solved most of
  %  the failed segmentations of subjects with fatty/thick skulls (ADNI,
  %  dll>0.1), whereas UKB subjects had dll<0.09. Adults (e.g. MRART,
  %  UKB) showed a log-likelihood gain for the children TPM of up to 0.07
  %  and children of 0.01-0.11. Further cases with failed segmentation are
  %  handled by the fallback in boney_segment_preprocessing.
  def.dllChildren  = 0.08; % minimal log-likelihood gain of the children TPM
  def.zoomChildren = 1.05; % minimal scaling (small head) for children
  def.dllBlaiotta  = 0.10; % minimal log-likelihood gain of the BlaiottaTPM

  sel.Ptpm = boney_segment_preprocessing('tpms');
  sel.Ptpm = sel.Ptpm(1:3); % the optional 4th TPM is not part of the selection
  V        = spm_vol(Pimg);

  % initial affine registration with the SPM TPM by the header origin or
  % the center of the field of view (see spm_preproc_run)
  tpm      = spm_load_priors8(sel.Ptpm{1});
  V1       = V; M = V1.mat; c = (V1.dim+1)/2; V1.mat(1:3,4) = -M(1:3,1:3)*c(:);
  evalc('[A1,ll1] = spm_maff8(V1,8,16,tpm,[],''subj'');'); A1 = A1*(V1.mat/M);
  evalc('[A2,ll2] = spm_maff8(V ,8,16,tpm,[],''subj'');');
  if ll1 > ll2, Ainit = A1; else, Ainit = A2; end

  % fast affine registration of all TPMs
  sel.ll = nan(1,numel(sel.Ptpm)); A = cell(1,numel(sel.Ptpm));
  for ti = 1:numel(sel.Ptpm)
    tpm = spm_load_priors8(sel.Ptpm{ti});
    evalc('[A{ti},sel.ll(ti)] = spm_maff8(V,6,16,tpm,Ainit,''subj'');');
  end
  sel.dll  = sel.ll - sel.ll(1);
  sel.zoom = abs(det(A{1}(1:3,1:3)))^(1/3);

  % selection
  if sel.dll(2) > def.dllChildren  &&  sel.zoom > def.zoomChildren
    ctpm = 2;
  elseif sel.dll(3) > def.dllBlaiotta  &&  sel.dll(3) > sel.dll(2)
    ctpm = 3;
  else
    ctpm = 1;
  end
  sel.ctpm = ctpm;

  % age-like parameter (only for humans, not a measure):
  % children have smaller heads (zoom>1) and fit better to the children
  % TPM, whereas older subjects fit better to the BlaiottaTPM (rough
  % linear model fitted on group ages of the development set, r=0.74,
  % mean absolute error ~10 years, children are clearly overestimated)
  sel.agelike = min(100, max(0, 78 - 175 * (sel.zoom - 1) - 375 * max(0,sel.dll(2)) + 30 * max(0,sel.dll(3))));
  sel.time    = etime(clock,stime);

  if verb
    tpmnames = {'SPM TPM','children TPM','BlaiottaTPM'};
    fprintf('  TPM selection: %s (agelike~%0.0f, zoom=%0.2f, dll=[%s], %0.0fs)\n', ...
      tpmnames{ctpm}, sel.agelike, sel.zoom, sprintf('%0.3f ',sel.dll(2:end)), sel.time);
  end
end
